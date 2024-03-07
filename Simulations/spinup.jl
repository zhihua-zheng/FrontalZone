### Setup dependencies
using Pkg; Pkg.instantiate()

#import NCDatasets as NCD

using Oceananigans
using Oceananigans.Units
using Oceananigans.Operators
using Oceananigans.TurbulenceClosures
using Oceananigans.BuoyancyModels: g_Earth
using Random, Printf

###########-------- SIMULATION PARAMETERS ----------------#############
const Lx = 1kilometers # east-west extent
const Ly = 1kilometers # north-south extent
const Lz = 140meters   # depth

const Nh_full = 512 # number of points in each of horizontal directions for full simulation
const Nz_full = 64  # number of points in the vertical direction for full simulation
const coarsen_factor_h = 1
const coarsen_factor_z = 1
Nx = Ny = Nh_full ÷ coarsen_factor_h
Nz = Nz_full ÷ coarsen_factor_z

const N₀² = 9e-8 # [s⁻²] mixed layer buoyancy frequency / stratification
const N₁² = 20*N₀² # [s⁻²] thermocline buoyancy frequency / stratification
const M² = 3e-8 # [s⁻²] horizontal buoyancy gradient
const f = 1e-4 # [s⁻¹] Coriolis frequency


###########-------- GRID SET UP ----------------#############
@info "Set up grid...."
const refinement = 1.25 # controls spacing near surface (higher means finer spaced)
const stretching = 8   # controls rate of stretching at bottom

# Normalized height ranging from 0 to 1
@inline h(k) = (k - 1) / Nz

# Linear near-surface generator
@inline ζ₀(k) = 1 + (h(k) - 1) / refinement

# Bottom-intensified stretching function
@inline Σ(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))

# Generating function
@inline z_faces(k) = Lz * (ζ₀(k) * Σ(k) - 1)

grid = RectilinearGrid(GPU(),
                       size = (Nx, Ny, Nz),
                       x = (-Lx/2, Lx/2),
                       y = (0, Ly),
                       z = z_faces,
                       topology = (Periodic, Periodic, Bounded))


###########-------- BOUNDARY CONDITIONS -----------------#############
@info "Set up boundary conditions...."
const Q₀ = 0.0   # [W m⁻²], surface heat flux (positive out of ocean)
const ρ₀ = 1026.0 # [kg m⁻³], average density at the surface of the world ocean
const cₚ = 3991.0  # [J K⁻¹ kg⁻¹], typical heat capacity for seawater
const αᵀ = 2e-4 # [K⁻¹], thermal expansion coefficient
const ν₀ = 1.0e-6 # [m² s⁻¹] molecular viscosity
const κ₀ = 1.5e-7 # [m² s⁻¹] molecular diffusivity
const hᵢ = 60 # [m] initial mixed layer depth

B₀ = g_Earth*αᵀ*Q₀ / (ρ₀*cₚ) # [m² s⁻³], surface buoyancy flux
const ∂v∂z0 = 0

u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0),
                                bottom = GradientBoundaryCondition(0))
v_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(-ν₀*∂v∂z0),
                                bottom = GradientBoundaryCondition(∂v∂z0))
w_bcs = FieldBoundaryConditions()
b_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(B₀),
                                bottom = GradientBoundaryCondition(N₁²))
eddy_νκ_bcs = FieldBoundaryConditions(top = ValueBoundaryCondition(0))


###########-------- TIME-INVARIANT BACKGROUND FIELDS -----------------#############
@info "Set up background fields...."
@inline B̅(x, y, z, t) = -M² * x
@inline V̅(x, y, z, t) = -M² / f * (z + hᵢ)

B_field = BackgroundField(B̅)
V_field = BackgroundField(V̅)


###########-------- INITIAL CONDITIONS -----------------#############
@info "Set up initial conditions...."
Random.seed!(45)
const noise = 1e-3
Ξ(z) = randn() * exp(4z/hᵢ)
uᵢ(x, y, z) = noise*Ξ(z)
vᵢ(x, y, z) = noise*Ξ(z)
wᵢ(x, y, z) = noise*Ξ(z)
bᵢ(x, y, z) = N₁²*(z + Lz) + (N₀² - N₁²)*max(z + hᵢ, 0)


###########-------- SPONGE LAYER -----------------#############
@info "Set up bottom sponge layer...."
# relax to initial states
const damping_rate = 1/60 #0.527*√N₁² # relax fields on a time-scale comparable to N₁, following Taylor & Ferrari 2010 
target_b = LinearTarget{:z}(intercept=N₁²*Lz, gradient=N₁²)
#@inline heaviside(X) = ifelse(X < 0, zero(X), one(X))
#@inline mask2nd(X) = heaviside(X) * X^2
#@inline function bottom_mask(x, y, z)
#    z₁ = -Lz; z₀ = z₁ + 20
#    return mask2nd((z₀ - z)/(z₀ - z₁))
#end
#@inline function bottom_mask(x, y, z)
#    z₁ = -Lz; z₀ = z₁ + 20; zₘ = (z₀ + z₁)/2
#    return (viside(z₀ - z) * (1 + tanh(3*(zₘ - z)/(zₘ - z₁)))/2
#end

#const sponge_z₁ = -Lz
#const sponge_z₀ = sponge_z₁ + 20
const sponge_σ = 6
#@inline bottom_mask(x, y, z) = heaviside(sponge_z₀ - z) * exp(-(z - sponge_z₁)^2 / (2*sponge_σ^2))
bottom_mask = GaussianMask{:z}(center=-Lz, width=sponge_σ)

uvw_sponge = Relaxation(rate=damping_rate, mask=bottom_mask, target=0)
b_sponge = Relaxation(rate=damping_rate, mask=bottom_mask, target=target_b)
sponge_forcing = (u=uvw_sponge, v=uvw_sponge, w=uvw_sponge, b=b_sponge)


###########-------- STARTING UP MODEL/ICs ---------------#############
@info "Define the model...."
#vitd = VerticallyImplicitTimeDiscretization()
model = NonhydrostaticModel(; grid,
                            coriolis = FPlane(f=f),
                            buoyancy = BuoyancyTracer(),
                            tracers = :b,
                            boundary_conditions = (b=b_bcs, u=u_bcs, v=v_bcs, w=w_bcs, νₑ=eddy_νκ_bcs, κₑ=(; b=eddy_νκ_bcs)),
                            forcing = sponge_forcing,
                            background_fields = (v=V_field, b=B_field),
                            advection = WENO(),
                            timestepper = :RungeKutta3,
                            closure = (ScalarDiffusivity(ν=ν₀, κ=κ₀), SmagorinskyLilly()))

set!(model, u=uᵢ, v=vᵢ, w=wᵢ, b=bᵢ)


###########-------- SIMULATION SET UP ---------------#############
@info "Define the simulation...."
Δx = minimum_xspacing(grid, Center(), Center(), Center())
Δy = minimum_yspacing(grid, Center(), Center(), Center())
Δz = minimum_zspacing(grid, Center(), Center(), Center())

const cfl_large = 0.9
Δt₀ = cfl_large * min(Δx, Δy, Δz) / 0.024

Tinertial = 2π/f
simulation = Simulation(model, Δt=Δt₀, stop_time=6*Tinertial, wall_time_limit=12hours)

wizard = TimeStepWizard(cfl=cfl_large, diffusive_cfl=cfl_large, min_change=0.05, max_change=1.5, max_Δt=5minutes)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(2))

wall_clock = Ref(time_ns())

function print_progress(sim)
    u, v, w = model.velocities
    progress = 100 * (time(sim) / sim.stop_time)
    elapsed = (time_ns() - wall_clock[]) / 1e9

    @printf("[%05.2f%%] i: %d, t: %s, wall time: %s, max(u): (%.1e, %.1e, %.1e) m/s, next Δt: %s\n",
            progress, iteration(sim), prettytime(sim), prettytime(elapsed),
            maximum(abs, u), maximum(abs, v), maximum(abs, w), prettytime(sim.Δt))

    wall_clock[] = time_ns()
    return nothing
end

simulation.callbacks[:print_progress] = Callback(print_progress, IterationInterval(100))


###########-------- CHECKPOINTER --------------#############
@info "Add checkpointer..."
ckpdir = "/glade/work/zhihuaz/Restart/FrontalZone/Front"
ispath(ckpdir) && rm(ckpdir, recursive=true, force=true)
mkdir(ckpdir)
simulation.output_writers[:checkpointer] = Checkpointer(model, schedule=TimeInterval(Tinertial), 
                                                        dir=ckpdir, prefix="checkpoint")


###########-------- RUN! --------------#############
run(`nvidia-smi`) # check how much memory used on a GPU run
@info "Run...."
run!(simulation)#, pickup="/glade/work/zhihuaz/Restart/FrontalZone/spinup/checkpoint_iteration17528.jld2")
@info "Simulation completed in " * prettytime(simulation.run_wall_time)
