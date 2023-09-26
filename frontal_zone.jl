### Setup dependencies
using Pkg; Pkg.instantiate()

using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures
using Printf

###########-------- SIMULATION PARAMETERS ----------------#############
const Lx = 1kilometers # east-west extent
const Ly = 1kilometers # north-south extent
const Lz = 140meters   # depth

const Nh_full = 512 # number of points in each of horizontal directions for full simulation
const Nz_full = 64  # number of points in the vertical direction for full simulation
const coarsen_factor = 4
Nx = Ny = Nh_full ÷ coarsen_factor
Nz = Nz_full ÷ coarsen_factor

const N₀² = 9e-8 # [s⁻²] mixed layer buoyancy frequency / stratification
const N₁² = 20*N₀² # [s⁻²] thermocline buoyancy frequency / stratification
const M² = 3e-8 # [s⁻²] horizontal buoyancy gradient
const f = 1e-4 # [s⁻¹] Coriolis frequency


###########-------- GRID SET UP ----------------#############
@info "Set up grid...."
const refinement = 1.2 # controls spacing near surface (higher means finer spaced)
const stretching = 12  # controls rate of stretching at bottom

# Normalized height ranging from 0 to 1
h(k) = (k - 1) / Nz

# Linear near-surface generator
ζ₀(k) = 1 + (h(k) - 1) / refinement

# Bottom-intensified stretching function
Σ(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))

# Generating function
z_faces(k) = Lz * (ζ₀(k) * Σ(k) - 1)

grid = RectilinearGrid(GPU(),
                       size = (Nx, Ny, Nz),
                       x = (-Lx/2, Lx/2),
                       y = (0, Ly),
                       z = z_faces,
                       topology = (Periodic, Periodic, Bounded))

###########-------- TIME-INVARIANT BACKGROUND FIELDS -----------------#############
@info "Set up background fields...."
parameters = (M2=M², f=f, H=Lz)
V(x, y, z, t, p) = p.M2 / p.f * (z + p.H)
B(x, y, z, t, p) = -p.M2 * x

V_field = BackgroundField(V, parameters=parameters)
B_field = BackgroundField(B, parameters=parameters)


###########-------- BOUNDARY CONDITIONS -----------------#############
@info "Set up boundary conditions...."
const Q₀ = 10.0   # W m⁻², surface heat flux (positive out of ocean)
const ρ₀ = 1026.0 # kg m⁻³, average density at the surface of the world ocean
const cₚ = 3991.0  # J K⁻¹ kg⁻¹, typical heat capacity for seawater
const g = 9.81 # m s⁻², gravitational acceleration
const αᵀ = 2e-4 # K⁻¹, thermal expansion coefficient
B₀ = g*αᵀ*Q₀ / (ρ₀*cₚ) # m² s⁻³, surface buoyancy flux
u_bcs = FieldBoundaryConditions(top = GradientBoundaryCondition(0),
                                bottom = GradientBoundaryCondition(0))
v_bcs = FieldBoundaryConditions(top = GradientBoundaryCondition(-M²/f),
                                bottom = GradientBoundaryCondition(-M²/f))
w_bcs = FieldBoundaryConditions()
b_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(B₀),
                                bottom = GradientBoundaryCondition(N₁²))


###########-------- INITIAL CONDITIONS -----------------#############
@info "Set up initial conditions...."
const noise = 1e-3
const hₘ = 60 # [m] initial mixed layer depth
uᵢ(x, y, z) = noise*randn()
vᵢ(x, y, z) = noise*randn()
wᵢ(x, y, z) = noise*randn()
bᵢ(x, y, z) = N₁²*(z + Lz) + (N₀² - N₁²)*max(z + hₘ, 0)


###########-------- SPONGE LAYER -----------------#############
@info "Set up bottom sponge layer...."
# relax to initial states
const damping_rate = 1/200 # relax fields on a 200 second time-scale
target_b = LinearTarget{:z}(intercept=N₁²*Lz, gradient=N₁²)
# bottom_mask = GaussianMask{:z}(center=-grid.Lz, width=20)
@inline heaviside(X) = ifelse(X < 0, zero(X), one(X))
@inline mask2nd(X) = heaviside(X) * X^2
@inline function bottom_mask(x, y, z)
    z₁ = -Lz; z₀ = z₁ + 20
    return mask2nd((z₀ - z)/(z₀ - z₁))
end

uvw_sponge = Relaxation(rate=damping_rate, mask=bottom_mask, target=0)
b_sponge = Relaxation(rate=damping_rate, mask=bottom_mask, target=target_b)
sponge_forcing = (u=uvw_sponge, v=uvw_sponge, w=uvw_sponge, b=b_sponge)


###########-------- STARTING UP MODEL/ICs ---------------#############
@info "Define the model...."
model = NonhydrostaticModel(; grid,
                            coriolis = FPlane(f=f),
                            buoyancy = BuoyancyTracer(),
                            tracers = :b,
                            boundary_conditions = (b=b_bcs, u=u_bcs, v=v_bcs, w=w_bcs),
                            forcing = sponge_forcing,
                            background_fields = (v=V_field, b=B_field),
                            advection = WENO(),
                            timestepper = :RungeKutta3,
                            closure = SmagorinskyLilly())

set!(model, u=uᵢ, v=vᵢ, w=wᵢ, b=bᵢ)


###########-------- SIMULATION SET UP ---------------#############
@info "Define the simulation...."
simulation = Simulation(model, Δt=20minutes, stop_time=20days)

wizard = TimeStepWizard(cfl=0.2, max_change=1.1, max_Δt=20minutes)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(20))

wall_clock = Ref(time_ns())

function print_progress(sim)
    u, v, w = model.velocities
    progress = 100 * (time(sim) / sim.stop_time)
    elapsed = (time_ns() - wall_clock[]) / 1e9

    @printf("[%05.2f%%] i: %d, t: %s, wall time: %s, max(u): (%6.3e, %6.3e, %6.3e) m/s, next Δt: %s\n",
            progress, iteration(sim), prettytime(sim), prettytime(elapsed),
            maximum(abs, u), maximum(abs, v), maximum(abs, w), prettytime(sim.Δt))

    wall_clock[] = time_ns()
    return nothing
end

simulation.callbacks[:print_progress] = Callback(print_progress, IterationInterval(100))


###########-------- DIAGNOSTICS --------------#############
@info "Add diagnostics..."
u, v, w = model.velocities
b = model.tracers.b
ζ = ∂x(v) - ∂y(u)
B = Average(b, dims=2)
U = Average(u, dims=2)
V = Average(v, dims=2)
W = Average(w, dims=2)
fields_slice = Dict("u" => u, "v" => v, "w" => w, "b" => b, "ζ" => ζ)
fields_zonnal_mean = Dict("B" => B, "U" => U, "V" => V, "W" => W)

filename = "frontal_zone"
save_fields_interval = 0.5day

slicers = (east = (grid.Nx, :, :),
           south = (:, 1, :),
           top = (:, :, grid.Nz))

for side in keys(slicers)
    indices = slicers[side]

    simulation.output_writers[side] = NetCDFOutputWriter(model, fields_slice;
                                                       filename = filename * "_$(side)_slice.nc",
                                                       schedule = TimeInterval(save_fields_interval),
                                                       overwrite_existing = true,
                                                       indices)
end

simulation.output_writers[:meridional] = NetCDFOutputWriter(model, fields_zonnal_mean;
                                                     filename = filename * "_meridional_mean.nc",
                                                     schedule = TimeInterval(save_fields_interval),
                                                     overwrite_existing = true)


###########-------- RUN! --------------#############
run(`nvidia-smi`) # check how much memory you're using on a GPU run
@info "Run...."
run!(simulation)
@info "Simulation completed in " * prettytime(simulation.run_wall_time)
