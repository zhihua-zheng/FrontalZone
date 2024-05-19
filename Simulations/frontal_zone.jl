### Setup dependencies
using Pkg; Pkg.instantiate()

using ArgParse
using Oceananigans
using Oceananigans.Units
using Oceananigans.Operators
using Oceananigans.TurbulenceClosures
using Oceananigans.BuoyancyModels: g_Earth
using Printf, Random, NCDatasets

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
           # default = "/glade/work/zhihuaz/Data/FrontalZone"   
            default = "/glade/derecho/scratch/zhihuaz/FrontalZone/Output"
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
groupname = ifelse(startswith(casename, 'f'), "Front", 
                   startswith(casename, 's') ? "ShortFront" : "NoFront")
pm = getproperty(SimParams(), Symbol(groupname))
pm = enrich_parameters(pm, casename)

stop_time = ifelse(args["nTf"] == nothing, pm.nTf, args["nTf"])*pm.Tf
ckpdir    = replace(outdir, "Output" => "Restart") * "/" * pm.ckpdir_affix


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
@inline V̅(x, y, z, t, p)  = -p.M² / p.f * (z + p.Lz/2)
@inline B̅(x, y, z, t, p)  = -p.M² * x

V_field = BackgroundField(V̅, parameters=background_params)
B_field = BackgroundField(B̅, parameters=background_params)
bkg_fields = ifelse(pm.use_background_Vg, (v=V_field, b=B_field), (; b=B_field))


###########-------- BOUNDARY CONDITIONS -----------------#############
@info "Set up boundary conditions...."
mom_flux_params = (; pm.σ_wind, pm.t₀, pm.tᵣ, pm.ν₀, pm.∂v∂z_cgeo, pm.τ₀ˣ, pm.τ₀ʸ, pm.ρ₀)
#@inline linear_ramp(t, p) = max((t - p.t₀), 0) / p.tᵣ
@inline cosine_ramp(t, p) = (1 - cos(π*max((t - p.t₀), 0) / p.tᵣ)) / 2
@inline ramp2const(t, p) = ifelse(t<(p.t₀+p.tᵣ), cosine_ramp(t, p), 1)
@inline oscillator(t, p) = sin(p.σ_wind * max((t - p.t₀), 0))
@inline temporal_wind(t, p) = ifelse(p.σ_wind!=0, oscillator(t, p), ramp2const(t, p))
@inline Fᵘ(x, y, t, p) = -p.τ₀ˣ/p.ρ₀*temporal_wind(t, p)
@inline Fᵛ(x, y, t, p) = -p.τ₀ʸ/p.ρ₀*temporal_wind(t, p) - p.ν₀*p.∂v∂z_cgeo

∂v∂z_bot = ifelse(pm.use_background_Vg, pm.∂v∂z_cgeo, -pm.M² / pm.f)
u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Fᵘ, parameters=mom_flux_params),
                                bottom = GradientBoundaryCondition(0))
v_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Fᵛ, parameters=mom_flux_params),
                                bottom = GradientBoundaryCondition(∂v∂z_bot))
w_bcs = FieldBoundaryConditions()
b_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(pm.B₀),
                                bottom = GradientBoundaryCondition(pm.N₁²))
#eddy_νκ_bcs = FieldBoundaryConditions(top = ValueBoundaryCondition(0))
eddy_νκ_bcs = FieldBoundaryConditions()


###########-------- SPONGE LAYER -----------------#############
@info "Set up bottom sponge layer...."
# relax to initial states
target_v = ifelse(pm.use_background_Vg, 0, LinearTarget{:z}(intercept=(-pm.M²/pm.f*pm.hᵢ),
                                                            gradient=(-pm.M²/pm.f)))
target_b = LinearTarget{:z}(intercept=pm.N₁²*pm.Lz, gradient=pm.N₁²)
#@inline heaviside(X) = ifelse(X < 0, zero(X), one(X))
#const sponge_z₁ = -Lz
#const sponge_z₀ = sponge_z₁ + 20
#@inline bottom_mask(x, y, z) = heaviside(sponge_z₀ - z) * exp(-(z - sponge_z₁)^2 / (2*sponge_σ^2))
bottom_mask = GaussianMask{:z}(center=-pm.Lz, width=pm.sponge_σ)

uw_sponge = Relaxation(rate=pm.damping_rate, mask=bottom_mask, target=0)
v_sponge  = Relaxation(rate=pm.damping_rate, mask=bottom_mask, target=target_v)
b_sponge  = Relaxation(rate=pm.damping_rate, mask=bottom_mask, target=target_b)
sponge_forcing = (u=uw_sponge, v=v_sponge, w=uw_sponge, b=b_sponge)


###########-------- STOKES DRIFT FORCING ---------------#############
@info "Set up Stokes drift forcing...."
Stokes_params = (; pm.Uˢ, pm.Dˢ, pm.θ₀, pm.t₀, pm.tᵣ)
@inline ∂z_uˢ(z, t, p) = p.Uˢ * exp(z / p.Dˢ) / p.Dˢ * cosd(p.θ₀)# * temporal_wind(t, p)
@inline ∂z_vˢ(z, t, p) = p.Uˢ * exp(z / p.Dˢ) / p.Dˢ * sind(p.θ₀)# * temporal_wind(t, p)

# layer-averaged Stokes drift
include("utils.jl")
#usla, vsla, dusdzla, dvsdzla = get_layer_averaged_stokes(grid)
fields_time_invariant = get_time_invariant_fields(grid)


###########-------- STARTING UP MODEL/ICs ---------------#############
@info "Define the model...."
model = NonhydrostaticModel(; grid,
                            coriolis = FPlane(f=pm.f),
                            buoyancy = BuoyancyTracer(),
                            tracers = (:b, :c),
                            stokes_drift = UniformStokesDrift(∂z_uˢ=∂z_uˢ, ∂z_vˢ=∂z_vˢ, parameters=Stokes_params),
                            boundary_conditions = (b=b_bcs, u=u_bcs, v=v_bcs, w=w_bcs, νₑ=eddy_νκ_bcs, κₑ=(; b=eddy_νκ_bcs)),
                            forcing = sponge_forcing,
                            background_fields = bkg_fields, 
                            advection = CenteredSecondOrder(),#WENO(),
                            timestepper = :RungeKutta3,
                            closure = SmagorinskyLilly())#ScalarDiffusivity(ν=pm.ν₀, κ=pm.κ₀),

if pm.pickup_checkpoint
    @info "Initialize from checkpoint file...."
    ckp_list  = split(read(`ls $ckpdir -1v`, String))
    ckp_fpath = ckpdir * "/" * ckp_list[7]
    set!(model, ckp_fpath)

    if pm.start_from_restratified
        @info "Initialize buoyancy from a restratified profile in frontal case..."
        bprof = sum(model.tracers.b.data.parent, dims=(1,2)) / (grid.Nx * grid.Ny)
        model.tracers.b.data.parent .-= bprof
        ds = NCDataset("/glade/work/zhihuaz/Data/FrontalZone/f11_Q000_W037_D000_St0_b-homo-restart.nc")
        bprof0 = ds["b"][:]
        copyto!(bprof, bprof0)
        model.tracers.b.data.parent .+= bprof
    end

    #file = jldopen(ckp_fpath, "r")
    #uᵢ = file["u/data"]
    #vᵢ = file["v/data"]
    #wᵢ = file["w/data"]
    #bᵢ = file["b/data"]
    #set!(model, u=uᵢ+usla, v=vᵢ+vsla, w=wᵢ, b=bᵢ, c=cᵢ)
else
    @info "Initialize from noise...."
    Random.seed!(45)
    Ξ(z) = randn() * exp(4z/pm.hᵢ)
    uᵢ(x, y, z) = pm.noise*Ξ(z)
    vᵢ(x, y, z) = ifelse(pm.use_background_Vg, pm.noise*Ξ(z), V̅(z)) 
    wᵢ(x, y, z) = pm.noise*Ξ(z)
    bᵢ(x, y, z) = pm.N₁²*(z + pm.Lz) + (pm.N₀² - pm.N₁²)*max(z + pm.hᵢ, 0)
    set!(model, u=uᵢ, v=vᵢ, w=wᵢ, b=bᵢ)
end

@info "Release tracer...."
cᵢ(x, y, z) = (1 - tanh(50*(z + pm.hᵢ) / pm.hᵢ)) / 2
set!(model, c=cᵢ)


###########-------- SIMULATION SET UP ---------------#############
@info "Define the simulation...."
Δx  = minimum_xspacing(grid, Center(), Center(), Center())
Δy  = minimum_yspacing(grid, Center(), Center(), Center())
Δz  = minimum_zspacing(grid, Center(), Center(), Center())
Δt₀ = pm.cfl * min(Δx, Δy, Δz) / max(pm.Vg, 0.02)
simulation = Simulation(model, Δt=Δt₀, stop_time=stop_time, wall_time_limit=12hours)

wizard = TimeStepWizard(cfl=pm.cfl, diffusive_cfl=pm.cfl, min_change=0.05, max_change=1.5, max_Δt=pm.max_Δt)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(2))

#wall_clock = Ref(time_ns())
#
#@inline function print_progress(sim)
#    u, v, w = model.velocities
#    progress = 100 * (time(sim) / sim.stop_time)
#    elapsed = (time_ns() - wall_clock[]) / 1e9
#
#    @printf("[%05.2f%%] i: %d, t: %s, wall time: %s, max(u): (%6.3e, %6.3e, %6.3e) m/s, next Δt: %s\n",
#            progress, iteration(sim), prettytime(sim), prettytime(elapsed),
#            maximum(abs, u), maximum(abs, v), maximum(abs, w), prettytime(sim.Δt))
#
#    wall_clock[] = time_ns()
#    return nothing
#end
#
#simulation.callbacks[:print_progress] = Callback(print_progress, IterationInterval(200))


###########-------- DIAGNOSTICS --------------#############
@info "Add diagnostics..."
include("diagnostics.jl")
fields_slice, fields_mean = get_output_tuple(model, fields_time_invariant["us"], fields_time_invariant["vs"],
                                             fields_time_invariant["Vbak"], pm; extra_outputs=pm.extra_outputs)

global_attributes = Dict("viscosity_mol" => pm.ν₀, "diffusivity_mol" => pm.κ₀,
                         "Uˢ" => pm.Uˢ, "Dˢ" => pm.Dˢ, "M²" => pm.M², "f" => pm.f)

#depth = 10:10:60
#depth_sym = (Symbol("xy", n) for n=depth)
#depth_idx = ((:, :, k) for k=[60, 54, 49, 43, 37, 32]) #XYSlice(grid, -depth)
#horizontal_slices = (; zip(depth_sym, depth_idx)...)
#vertical_slices   = (yz = (grid.Nx, :, :),
#                     xz = (:, 1, :))
#selected_slices   = merge(vertical_slices, horizontal_slices)
#slicers = pm.full_fields ? (full = (:,:,:),) : selected_slices
slicers = (full = (:,:,:),)

for side in keys(slicers)
    indices = slicers[side]
    simulation.output_writers[side] = NetCDFOutputWriter(model, fields_slice;
                                                       filename = casename * "_$(side).nc",
                                                       dir = outdir,
                                                       schedule = TimeInterval(pm.out_interval_slice),
                                                       global_attributes = global_attributes,
                                                       overwrite_existing = true,
                                                       indices)
end

ow = simulation.output_writers[:averages] = NetCDFOutputWriter(model, fields_mean;
                                                       filename = casename * "_averages.nc",
                                                       dir = outdir,
                                                       schedule = TimeInterval(pm.out_interval_mean),
                                                       global_attributes = global_attributes,
                                                       overwrite_existing = true)

write_time_invariant_fields!(model, ow, fields_time_invariant)


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
