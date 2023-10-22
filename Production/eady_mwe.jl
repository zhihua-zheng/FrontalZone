using Random, Printf, Oceananigans, Oceananigans.Units, Oceananigans.TurbulenceClosures

Random.seed!(43)

h(k) = (k - 1) / Nz
ζ₀(k) = 1 + (h(k) - 1) / refinement
Σ(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))
z_faces(k) = Lz * (ζ₀(k) * Σ(k) - 1)

Nx, Ny, Nz = 512, 512, 64
Lx, Ly, Lz = 1kilometer, 1kilometer, 140.0

refinement = 1.8
stretching = 3

arch = GPU()
grid = RectilinearGrid(arch; size = (Nx, Ny, Nz), x = (0, Lx), y = (0, Ly), z = z_faces)

coriolis = FPlane(f = 1e-4) # [s⁻¹]
background_state_parameters = (M = 1e-4,       # s⁻¹, geostrophic shear
                               f = coriolis.f, # s⁻¹, Coriolis parameter
                               N = 1e-4,       # s⁻¹, buoyancy frequency
                               H = grid.Lz )
B(x, y, z, t, p) = p.M ^ 2 * x + p.N ^ 2 * (z + p.H)
V(x, y, z, t, p) = p.M ^ 2 / p.f * (z + p.H)

V_field = BackgroundField(V, parameters = background_state_parameters)
B_field = BackgroundField(B, parameters = background_state_parameters)

νᵥ = κᵥ = 1e-4 # [m² s⁻¹]
closure = ScalarDiffusivity(ν = νᵥ, κ = κᵥ)
#closure = SmagorinskyLilly() 

model = NonhydrostaticModel(; grid,
                            advection=WENO(grid),
                            timestepper = :RungeKutta3,
                            coriolis,
                            tracers = :b,
                            buoyancy = BuoyancyTracer(),
                            background_fields = (b = B_field, v = V_field),
                            closure)

model.clock.time = 50days

Ξ(z) = randn() * z / grid.Lz * (z / grid.Lz + 1)

Ũ = 1e-3
uᵢ(x, y, z) = Ũ * Ξ(z)
vᵢ(x, y, z) = Ũ * Ξ(z)

set!(model, u=uᵢ, v=vᵢ)

Δx = minimum_xspacing(grid, Center(), Center(), Center())
Δy = minimum_yspacing(grid, Center(), Center(), Center())
Δz = minimum_zspacing(grid, Center(), Center(), Center())

#const cfl_small = 0.01
const cfl_large = 0.75
#cfl(t) = cfl_small + (cfl_large - cfl_small) * (1 + tanh( (t - 52days) / 0.5days )) / 2
#max(min((t - 50days)/1days - 0.15, cfl_large), cfl_small)
Δt₀ = cfl_large * min(Δx, Δy, Δz) / V(0, 0, 0, 0, background_state_parameters)

simulation = Simulation(model, Δt = Δt₀, stop_time = 58days)

# Adapt the time step while keeping the CFL number fixed.
wizard = TimeStepWizard(cfl = cfl_large, diffusive_cfl = cfl_large, max_Δt = 3minutes, min_change = 0.1, max_change = 1.5)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(1))

#progress(sim) = @printf("i: % 6d, sim time: % 10s, wall time: % 10s, Δt: % 10s\n",
#                        sim.model.clock.iteration,
#                        prettytime(sim.model.clock.time),
#                        prettytime(sim.run_wall_time),
#                        prettytime(sim.Δt))
function progress(sim)
    sim.model.clock.time > 54days && (wizard.max_Δt = 5minutes)
    #wizard.cfl = cfl(model.clock.time)
    @printf("i: % 6d, sim time: % 10s, wall time: % 10s, Δt: % 10s, advective CFL: %.2e, diffusive CFL: %.2e\n",
            sim.model.clock.iteration,
            prettytime(sim.model.clock.time),
            prettytime(sim.run_wall_time),
            prettytime(sim.Δt),
            AdvectiveCFL(sim.Δt)(sim.model),
            DiffusiveCFL(sim.Δt)(sim.model))
end
simulation.callbacks[:progress] = Callback(progress, IterationInterval(1))

u, v, w = model.velocities
b = model.tracers.b
fields_slice = Dict("u" => u, "v" => v, "w" => w, "b" => b)
filename = "eady_mwe.nc"
data_dir = "/glade/work/zhihuaz/Data/FrontalZone"

simulation.output_writers[:fields] = NetCDFOutputWriter(model, fields_slice;
                                                     filename = filename,
                                                     dir = data_dir,
                                                     schedule = TimeInterval(6hours),
                                                     overwrite_existing = true)

run(`nvidia-smi`)
run!(simulation)
