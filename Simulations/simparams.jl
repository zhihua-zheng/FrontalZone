using Parameters

@with_kw struct SimParams
    commons = (;
               Lz = 140meters, # depth
               Nz_full = 64,   # number of points in the z direction for full simulation

               N₀² = 9e-8,   # [s⁻²] mixed layer buoyancy frequency (stratification)
               f   = 1e-4,   # [s⁻¹] Coriolis frequency
               hᵢ  = 60,     # [m] initial mixed layer depth
               ρₐ  = 1.225,  # [kg m⁻³] average density of air at the surface
               ρ₀  = 1026.0, # [kg m⁻³] average density of seawater at the surface
               cₚ  = 3991.0, # [J K⁻¹ kg⁻¹] typical heat capacity for seawater
               αᵀ  = 2e-4,   # [K⁻¹], thermal expansion coefficient
               ν₀  = 1.0e-6, # [m² s⁻¹] molecular viscosity
               κ₀  = 1.5e-7, # [m² s⁻¹] molecular diffusivity

               noise = 1e-3, # initial noise amplitude
               out_interval = 1hour, # how often to write output

               z_refinement = 1.25, # controls spacing near surface (higher means finer spaced)
               z_stretching = 8,    # controls rate of stretching at bottom
               sponge_σ     = 6,    # [m] sponge layer Gaussian mask width
               damping_rate = 1/60, # [s⁻¹] relax fields on a time-scale comparable to N₁, following Taylor & Ferrari 2010
 
               # whether use compensating stress for perturbation along-front velocity in unforced conditions with background velocity
               stress_with_bgV_unforced = false,
               extra_outputs = false,
              )

    Front   = (; commons..., M² = 3e-8, # [s⁻²] horizontal buoyancy gradient
               Lx = 1kilometers, # east-west extent
               Ly = 1kilometers, # north-south extent
               Nx_full = 512, # number of points in the x direction for full simulation
               Ny_full = 512, # number of points in the y direction for full simulation
               nTf = 12,      # simulation length in unit of inertial period
               cfl = 0.9,     # Courant-Friedrichs–Lewy number
               )

    ShortFront = (; commons..., M² = 3e-8,
               Lx = 1kilometers,
               Ly = 250meters,
               Nx_full = 512,
               Ny_full = 128,
               nTf = 12,
               cfl = 0.9,
               )

    NoFront = (; commons..., M² = 0,
               Lx = 1kilometers,
               Ly = 1kilometers,
               Nx_full = 512,
               Ny_full = 512,
               nTf = 7.2,
               cfl = 0.65,
               )
end


function decode_casename(casename)
    fres, heat_flux, wind_stress, wind_direction, Stokes_flag = split(casename, "_")
    coarsen_factor_h = parse(Int64, fres[2])
    coarsen_factor_z = parse(Int64, fres[3])
    Q₀ = parse(Float64, lstrip(heat_flux,      'Q'))
    τ₀ = parse(Float64, lstrip(wind_stress,    'W')) / 1e3
    if startswith(wind_direction, 'D')
        θ₀ = parse(Float64, lstrip(wind_direction, 'D'))
        wind_oscillation = false
    elseif startswith(wind_direction, 'O')
        θ₀ = parse(Float64, lstrip(wind_direction, 'O'))
        wind_oscillation = true
    else
        throw(DomainError(wind_direction, "wind direction identifier is bad"))
    end
    use_Stokes = parse(Bool, Stokes_flag[end])
    return coarsen_factor_h, coarsen_factor_z, Q₀, τ₀, θ₀, wind_oscillation, use_Stokes
end


function enrich_parameters(params, casename)
    coarsen_factor_h, coarsen_factor_z, Q₀, τ₀, θ₀, wind_oscillation, use_Stokes = decode_casename(casename)
    σ_wind = ifelse(wind_oscillation, params.f, 0)
    save_checkpoint   = ifelse((Q₀ + τ₀)==0, true, false)
    pickup_checkpoint = ifelse(startswith(casename, 'n'), false, !save_checkpoint)

    # vertical gradient of along-front velocity at the boundary
    ∂v∂z_uf = ifelse(pm.stress_with_bgV_unforced, pm.M²/pm.f, 0)
    N₁² = 20*params.N₀² # [s⁻²] thermocline buoyancy frequency (stratification)
    Tf  = 2π/params.f   # [s] inertial period 
    t₀  = ifelse(startswith(casename, 'n'), 0, 5*Tf) # [s] when to introduce wind stres
    tᵣ  = 2*Tf # [s] length of the linear ramp of wind forcing
    ckp_interval = 1*Tf # [s] how often to checpoint

    Nx = params.Nx_full ÷ coarsen_factor_h
    Ny = params.Ny_full ÷ coarsen_factor_h
    Nz = params.Nz_full ÷ coarsen_factor_z

    τ₀ˣ   = τ₀ * cosd(θ₀)
    τ₀ʸ   = τ₀ * sind(θ₀)
    Cd    = 1.2e-3 # neutral drag coefficient (for U10 <= 11 m s⁻¹) from Large & Pond 1981
    U₁₀   = √(τ₀ / params.ρₐ / Cd) # [m s⁻¹] surface wind speed corresponding to τ₀
    Uˢ    = ifelse(use_Stokes, 0.0155 * U₁₀, 0) # [m s⁻¹] surface Stokes drift velocity
    Dˢ    = ifelse(use_Stokes, 0.14 * U₁₀^2 / g_Earth, 1) # [m] vertical (e-folding) scale of the Stokes drift (avoid blowing up)
    B₀    = g_Earth * params.αᵀ * Q₀ / (params.ρ₀ * params.cₚ) # [m² s⁻³] surface buoyancy flux
    ustar = √(τ₀ / params.ρ₀) # [m s⁻¹] friction velocity
    wstar = ∛(B₀ * params.hᵢ) # [m s⁻¹] convective velocity

    Vg   = params.M² / params.f * params.Lz # [m s⁻¹] geostrophic velocity scale
    RiB₀ = params.N₀² * params.f^2 / max(params.M²^2, 1e-24) # mixed layer balanced Richardson number
    RiB₁ = N₁² * params.f^2 / max(params.M²^2, 1e-24)        # thermocline ...

    extra_params = Base.@locals()
    delete!(extra_params, :params)
    delete!(extra_params, :casename)
    return merge(params, NamedTuple(extra_params))
end
