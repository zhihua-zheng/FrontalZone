### Setup dependencies
#using Pkg; Pkg.instantiate()

#import NCDatasets as NCD
using ArgParse
using Oceananigans
using Oceananigans.Units
using Oceananigans.Operators
using Oceananigans.TurbulenceClosures
using Oceananigans.BuoyancyModels: g_Earth
using Printf, Random

###########-------- COMMAND LINE ARGUMENTS ----------------#############
@info "Parse command line arguments..."
# Returns a dictionary of command line arguments
function parse_command_line_arguments()
    settings = ArgParseSettings()
    @add_arg_table! settings begin
        "casename"
            help = "Name of simulation case"
            required = true
            arg_type = String

        "--nTf"
            help = "Number of inertial periods the simulation runs"
            arg_type = Float64

        "--outdir"
            help = "Path of directory to save outputs under"
            default = "/glade/work/zhihuaz/Data/FrontalZone"   
            arg_type = String
    end
    return parse_args(settings)
end

args = parse_command_line_arguments()
for (arg,val) in args
    @info "    $arg => $val"
end

casename = args["casename"]
outdir   = args["outdir"]


###########-------- SIMULATION PARAMETERS ----------------#############
@info "Load in simulation parameters..."
include("simparams.jl")
groupname = ifelse(startswith(casename, 'f'), "Front", "NoFront")
pm = getproperty(SimParams(), Symbol(groupname))
pm = enrich_parameters(pm, casename)

stop_time = ifelse(args["nTf"]==nothing, pm.nTf, args["nTf"])*pm.Tf
ckpdir    = replace(outdir, "Data" => "Restart") * "/" * groupname


###########-------- GRID SET UP ----------------#############
@info "Set up grid...."
# Normalized height ranging from 0 to 1
@inline h(k) = (k - 1) / pm.Nz

# Linear near-surface generator
@inline ζ₀(k) = 1 + (h(k) - 1) / pm.z_refinement

# Bottom-intensified stretching function
@inline Σ(k) = (1 - exp(-pm.z_stretching * h(k))) / (1 - exp(-pm.z_stretching))

# Generating function
@inline z_faces(k) = pm.Lz * (ζ₀(k) * Σ(k) - 1)

grid = RectilinearGrid(GPU(),
                       size = (pm.Nx, pm.Ny, pm.Nz),
                       x = (-pm.Lx/2, pm.Lx/2),
                       y = (0, pm.Ly),
                       z = z_faces,
                       topology = (Periodic, Periodic, Bounded))


###########-------- TIME-INVARIANT BACKGROUND FIELDS -----------------#############
@info "Set up background fields...."
background_params = (; pm.M², pm.f, pm.Lz)
V̅(x, y, z, t, p)  = -p.M² / p.f * (z + p.Lz)
B̅(x, y, z, t, p)  = -p.M² * x

V_field = BackgroundField(V̅, parameters=background_params)
B_field = BackgroundField(B̅, parameters=background_params)


###########-------- BOUNDARY CONDITIONS -----------------#############
@info "Set up boundary conditions...."
mom_flux_params = (; pm.t₀, pm.tᵣ, pm.ν₀, pm.∂v∂z_uf, pm.τ₀ˣ, pm.τ₀ʸ, pm.ρ₀)
@inline linear_ramp(t, p) = min( max( (t - p.t₀)/p.tᵣ, 0 ), 1 )
@inline Fᵘ(x, y, t, p) = -linear_ramp(t, p) * p.τ₀ˣ/p.ρ₀
@inline Fᵛ(x, y, t, p) = -linear_ramp(t, p) * p.τ₀ʸ/p.ρ₀ - p.ν₀*p.∂v∂z_uf

u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Fᵘ, parameters=mom_flux_params),
                                bottom = GradientBoundaryCondition(0))
v_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Fᵛ, parameters=mom_flux_params),
                                bottom = GradientBoundaryCondition(pm.∂v∂z_uf))
w_bcs = FieldBoundaryConditions()
b_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(pm.B₀),
                                bottom = GradientBoundaryCondition(pm.N₁²))
eddy_νκ_bcs = FieldBoundaryConditions(top = ValueBoundaryCondition(0))


###########-------- SPONGE LAYER -----------------#############
@info "Set up bottom sponge layer...."
# relax to initial states
target_b = LinearTarget{:z}(intercept=pm.N₁²*pm.Lz, gradient=pm.N₁²)
#@inline heaviside(X) = ifelse(X < 0, zero(X), one(X))
#const sponge_z₁ = -Lz
#const sponge_z₀ = sponge_z₁ + 20
#@inline bottom_mask(x, y, z) = heaviside(sponge_z₀ - z) * exp(-(z - sponge_z₁)^2 / (2*sponge_σ^2))
bottom_mask = GaussianMask{:z}(center=-pm.Lz, width=pm.sponge_σ)

uvw_sponge = Relaxation(rate=pm.damping_rate, mask=bottom_mask, target=0)
b_sponge   = Relaxation(rate=pm.damping_rate, mask=bottom_mask, target=target_b)
sponge_forcing = (u=uvw_sponge, v=uvw_sponge, w=uvw_sponge, b=b_sponge)


###########-------- STOKES DRIFT FORCING ---------------#############
@info "Set up Stokes drift forcing...."
Stokes_params = (; pm.Uˢ, pm.Dˢ)
#uˢ(z, p) = p.Uˢ * exp(z / p.Dˢ)
∂z_uˢ(z, t, p) = p.Uˢ * exp(z / p.Dˢ) / p.Dˢ


###########-------- STARTING UP MODEL/ICs ---------------#############
@info "Define the model...."
model = NonhydrostaticModel(; grid,
                            coriolis = FPlane(f=pm.f),
                            buoyancy = BuoyancyTracer(),
                            tracers = :b,
                            stokes_drift = UniformStokesDrift(∂z_uˢ=∂z_uˢ, parameters=Stokes_params),
                            boundary_conditions = (b=b_bcs, u=u_bcs, v=v_bcs, w=w_bcs, νₑ=eddy_νκ_bcs, κₑ=(; b=eddy_νκ_bcs)),
                            forcing = sponge_forcing,
                            background_fields = (v=V_field, b=B_field),
                            advection = WENO(),
                            timestepper = :RungeKutta3,
                            closure = (ScalarDiffusivity(ν=pm.ν₀, κ=pm.κ₀), SmagorinskyLilly()))

if pm.pickup_checkpoint
    @info "Initialize from checkpoint file...."
    #ckp_list = split(read(`ls $ckpdir -1v`, String))
    #irestart = Int64(pm.t₀ / pm.ckp_interval) + 1
    ckp_file = ckpdir * "/" * "checkpoint_iteration7426.jld2"#ckp_list[irestart] 
    set!(model, ckp_file)
else
    @info "Initialize from noise...."
    Random.seed!(45)
    Ξ(z) = randn() * exp(4z/pm.hᵢ)
    uᵢ(x, y, z) = pm.noise*Ξ(z)
    vᵢ(x, y, z) = pm.noise*Ξ(z)
    wᵢ(x, y, z) = pm.noise*Ξ(z)
    bᵢ(x, y, z) = pm.N₁²*(z + pm.Lz) + (pm.N₀² - pm.N₁²)*max(z + pm.hᵢ, 0)
    set!(model, u=uᵢ, v=vᵢ, w=wᵢ, b=bᵢ)
end


###########-------- SIMULATION SET UP ---------------#############
@info "Define the simulation...."
Δx = minimum_xspacing(grid, Center(), Center(), Center())
Δy = minimum_yspacing(grid, Center(), Center(), Center())
Δz = minimum_zspacing(grid, Center(), Center(), Center())
Δt₀ = pm.cfl * min(Δx, Δy, Δz) / max(pm.Vg, 0.02)
simulation = Simulation(model, Δt=Δt₀, stop_time=stop_time, wall_time_limit=12hours)

wizard = TimeStepWizard(cfl=pm.cfl, diffusive_cfl=pm.cfl, min_change=0.05, max_change=1.5, max_Δt=5minutes)
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

simulation.callbacks[:print_progress] = Callback(print_progress, IterationInterval(1000))


###########-------- DIAGNOSTICS --------------#############
@info "Add diagnostics..."
include("diagnostics.jl")
fields_slice, fields_mean = get_output_tuple(model)

global_attributes = Dict("viscosity_mol" => pm.ν₀, "diffusivity_mol" => pm.κ₀, "M²" => pm.M², "f" => pm.f)
slicers = (east = (grid.Nx, :, :),
           south = (:, 1, :),
           top = (:, :, grid.Nz-5))

for side in keys(slicers)
    indices = slicers[side]
    simulation.output_writers[side] = NetCDFOutputWriter(model, fields_slice;
                                                       filename = casename * "_$(side)_slice.nc",
                                                       dir = outdir,
                                                       schedule = TimeInterval(pm.out_interval),
                                                       global_attributes = global_attributes,
                                                       overwrite_existing = true,
                                                       indices)
end

simulation.output_writers[:averages] = NetCDFOutputWriter(model, fields_mean;
                                                     filename = casename * "_averages.nc",
                                                     dir = outdir,
                                                     schedule = TimeInterval(pm.out_interval), 
                                                     global_attributes = global_attributes,
                                                     overwrite_existing = true)


###########-------- CHECKPOINTER --------------#############
if pm.save_checkpoint
    @info "Add checkpointer..."
    ispath(ckpdir) && rm(ckpdir, recursive=true, force=true) 
    mkdir(ckpdir)
    simulation.output_writers[:checkpointer] = Checkpointer(model, schedule=TimeInterval(pm.ckp_interval),
                                                            dir=ckpdir, prefix="checkpoint")
end


###########-------- RUN! --------------#############
run(`nvidia-smi`) # check how much memory used on a GPU run
@info "Run...."
run!(simulation)
@info "Simulation completed in " * prettytime(simulation.run_wall_time)
