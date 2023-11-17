### Setup dependencies
#using Pkg; Pkg.instantiate()

#import NCDatasets as NCD

using Oceananigans
using Oceananigans.Units
using Oceananigans.Operators
using Oceananigans.TurbulenceClosures
using Oceananigans.BuoyancyModels: g_Earth
using Printf

###########-------- SIMULATION PARAMETERS ----------------#############
casename = "r11-Q000-W022-D00-St0"
use_Stokes = false
save_checkpoint = false

# `noVflux_total` decides if the simulation applies no flux top BC for total velocity in unforced conditions
# If `noVflux_total = false`, no flux top BC is applied for perturbation velocity,
# implying a stress to maitain the background flow against dissipation.
noVflux_total = false

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

###########-------- TIME-INVARIANT BACKGROUND FIELDS -----------------#############
@info "Set up background fields...."
parameters = (M2=M², f=f, H=Lz)
V̅(x, y, z, t, p) = -p.M2 / p.f * (z + p.H)
B̅(x, y, z, t, p) = -p.M2 * x

V_field = BackgroundField(V̅, parameters=parameters)
B_field = BackgroundField(B̅, parameters=parameters)


###########-------- BOUNDARY CONDITIONS -----------------#############
@info "Set up boundary conditions...."
const τ₀ = 0.022 # [N m⁻²], surface wind stress
const θ₀ = 0.0 # [degree], surface wind direction, relative to x-axis, positive counter-clockwise
const Q₀ = 0.0 # [W m⁻²], surface heat flux (positive out of ocean)
const ρ₀ = 1026.0 # [kg m⁻³], average density at the surface of the world ocean
const cₚ = 3991.0 # [J K⁻¹ kg⁻¹], typical heat capacity for seawater
const αᵀ = 2e-4 # [K⁻¹], thermal expansion coefficient
const ν₀ = 1.0e-6 # [m² s⁻¹] molecular viscosity
const κ₀ = 1.5e-7 # [m² s⁻¹] molecular diffusivity
#const hᵢ = 60 # [m] initial mixed layer depth

B₀ = g_Earth*αᵀ*Q₀ / (ρ₀*cₚ) # [m² s⁻³], surface buoyancy flux
# This sets the surface gradient of along-front perturbation velocity v in unforced conditions
if noVflux_total
    ∂v∂z0 = M²/f
else
    ∂v∂z0 = 0
end
@inline linear_ramp(t, t₀, tᵣ) = min((t - t₀)/tᵣ, 1)
@inline Fᵘ(x, y, t, p) = -linear_ramp(t, p.t₀, p.tᵣ) * τ₀/ρ₀ * cosd(θ₀)
@inline Fᵛ(x, y, t, p) = -linear_ramp(t, p.t₀, p.tᵣ) * τ₀/ρ₀ * sind(θ₀) - ν₀*p.∂v∂z0
mom_flux_params = (t₀=3.5days, # time to introduce wind stress 
                   tᵣ=1.5days, # length of linear ramp
                   ∂v∂z0=∂v∂z0)

u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Fᵘ, parameters=mom_flux_params),
                                bottom = GradientBoundaryCondition(0))
v_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Fᵛ, parameters=mom_flux_params),
                                bottom = GradientBoundaryCondition(∂v∂z0))
w_bcs = FieldBoundaryConditions()
b_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(B₀),
                                bottom = GradientBoundaryCondition(N₁²))
eddy_νκ_bcs = FieldBoundaryConditions(top = ValueBoundaryCondition(0))


###########-------- SPONGE LAYER -----------------#############
@info "Set up bottom sponge layer...."
# relax to initial states
const damping_rate = 1/60 #0.527*√N₁² # relax fields on a time-scale comparable to N₁, following Taylor & Ferrari 2010
target_b = LinearTarget{:z}(intercept=N₁²*Lz, gradient=N₁²)
#@inline heaviside(X) = ifelse(X < 0, zero(X), one(X))
#const sponge_z₁ = -Lz
#const sponge_z₀ = sponge_z₁ + 20
const sponge_σ = 6
#@inline bottom_mask(x, y, z) = heaviside(sponge_z₀ - z) * exp(-(z - sponge_z₁)^2 / (2*sponge_σ^2))
bottom_mask = GaussianMask{:z}(center=-Lz, width=sponge_σ)

uvw_sponge = Relaxation(rate=damping_rate, mask=bottom_mask, target=0)
b_sponge = Relaxation(rate=damping_rate, mask=bottom_mask, target=target_b)
sponge_forcing = (u=uvw_sponge, v=uvw_sponge, w=uvw_sponge, b=b_sponge)


###########-------- STOKES DRIFT FORCING ---------------#############
@info "Set up Stokes drift forcing...."
#amplitude  = 0.8 # [m]
#wavelength = 60  # [m]
#wavenumber = 2π / wavelength # [m⁻¹]
#frequency  = sqrt(g_Earth * wavenumber) # [s⁻¹]

# The vertical scale over which the Stokes drift of a monochromatic surface wave
# decays away from the surface is `1/2wavenumber`, or
#const vertical_scale = wavelength / 4π

# Stokes drift velocity at the surface
#const Uˢ = amplitude^2 * wavenumber * frequency # m s⁻¹

if use_Stokes
    uˢ(z) = Uˢ * exp(z / vertical_scale)
    ∂z_uˢ(z, t) = 1 / vertical_scale * Uˢ * exp(z / vertical_scale)
else
    uˢ(z) = 0
    ∂z_uˢ(z, t) = 0
end


###########-------- STARTING UP MODEL/ICs ---------------#############
@info "Define the model...."
model = NonhydrostaticModel(; grid,
                            coriolis = FPlane(f=f),
                            buoyancy = BuoyancyTracer(),
                            tracers = :b,
                            stokes_drift = UniformStokesDrift(∂z_uˢ=∂z_uˢ),
                            boundary_conditions = (b=b_bcs, u=u_bcs, v=v_bcs, w=w_bcs, νₑ=eddy_νκ_bcs, κₑ=(; b=eddy_νκ_bcs)),
                            forcing = sponge_forcing,
                            background_fields = (v=V_field, b=B_field),
                            advection = WENO(),
                            timestepper = :RungeKutta3,
                            closure = (ScalarDiffusivity(ν=ν₀, κ=κ₀), SmagorinskyLilly()))

ckp_file = "/glade/work/zhihuaz/Restart/FrontalZone/spinup/checkpoint_iteration7867.jld2" 
set!(model, ckp_file)


###########-------- SIMULATION SET UP ---------------#############
@info "Define the simulation...."
Δx = minimum_xspacing(grid, Center(), Center(), Center())
Δy = minimum_yspacing(grid, Center(), Center(), Center())
Δz = minimum_zspacing(grid, Center(), Center(), Center())

const cfl_large = 0.9
Δt₀ = cfl_large * min(Δx, Δy, Δz) / abs(V̅(0, 0, 0, 0, parameters))

simulation = Simulation(model, Δt=Δt₀, stop_time=8.5days, wall_time_limit=12hours)

wizard = TimeStepWizard(cfl=cfl_large, diffusive_cfl=cfl_large, min_change=0.05, max_change=1.5, max_Δt=5minutes)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(2))

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

simulation.callbacks[:print_progress] = Callback(print_progress, IterationInterval(500))


###########-------- DIAGNOSTICS --------------#############
@info "Add diagnostics..."
include("diagnostics.jl")
fields_slice, fields_mean, fields_mean_extra = get_output_tuple(model)

global_attributes = Dict("viscosity_mol" => ν₀, "diffusivity_mol" => κ₀, "M²" => M², "f" => f)
data_dir = "/glade/work/zhihuaz/Data/FrontalZone"
save_fields_interval = 1hour

#zC = znodes(grid, Nothing, Nothing, Center())
#idx_z_slice = findfirst(n -> n === min(zC), zC)
slicers = (east = (grid.Nx, :, :),
           south = (:, 1, :),
           top = (:, :, grid.Nz-5))

for side in keys(slicers)
    indices = slicers[side]
    simulation.output_writers[side] = NetCDFOutputWriter(model, fields_slice;
                                                       filename = casename * "_$(side)_slice.nc",
                                                       dir = data_dir,
                                                       schedule = TimeInterval(save_fields_interval),
                                                       global_attributes = global_attributes,
                                                       overwrite_existing = true,
                                                       indices)
end

simulation.output_writers[:averages] = NetCDFOutputWriter(model, fields_mean;
                                                     filename = casename * "_averages.nc",
                                                     dir = data_dir,
                                                     schedule = TimeInterval(save_fields_interval), 
                                                     global_attributes = global_attributes,
                                                     overwrite_existing = true)


simulation.output_writers[:averages_extra] = NetCDFOutputWriter(model, fields_mean_extra;
                                                     filename = casename * "_averages_extra.nc",
                                                     dir = data_dir,
                                                     schedule = TimeInterval(1hour), 
                                                     global_attributes = global_attributes,
                                                     overwrite_existing = true)

###########-------- RUN! --------------#############
run(`nvidia-smi`) # check how much memory used on a GPU run
@info "Run...."
run!(simulation)
@info "Simulation completed in " * prettytime(simulation.run_wall_time)
