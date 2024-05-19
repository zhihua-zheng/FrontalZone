using Parameters

@with_kw struct SimParams
    commons = (;
               f  = 1e-4,   # [s⁻¹] Coriolis frequency
               hᵢ = 60,     # [m] initial mixed layer depth
               ρₐ = 1.225,  # [kg m⁻³] average density of air at the surface
               ρ₀ = 1026.0, # [kg m⁻³] average density of seawater at the surface
               cₚ = 3991.0, # [J K⁻¹ kg⁻¹] typical heat capacity for seawater
               αᵀ = 2e-4,   # [K⁻¹], thermal expansion coefficient
               ν₀ = 1.0e-6, # [m² s⁻¹] molecular viscosity
               κ₀ = 1.5e-7, # [m² s⁻¹] molecular diffusivity

               RiB₀ = 1,  # mixed layer balanced Richardson number
               RiB₁ = 20, # thermocline balanced Richardson number

               noise = 1e-3, # initial noise amplitude
               out_interval_slice = 15minutes, # how often to save slice
               out_interval_mean  = 0.5minutes, # how often to save mean

               z_refinement = 3.7,  # controls spacing near surface (higher means finer spaced)
               z_stretching = 4,  # controls rate of stretching at bottom
               damping_rate = 1/60, # [s⁻¹] relax fields on a time-scale comparable to N₁, following Taylor & Ferrari 2010
 
               # whether to apply additional forcing to cancel the geostrophic stress associated with the background velocity
               counter_geo_stress = false,
               use_background_Vg  = true,
               extra_outputs      = false,
               save_checkpoint    = false,
               full_fields        = true,
              )

    Front   = (; commons...,
               Lx = 1kilometers, # east-west extent
               Ly = 1kilometers, # north-south extent
               Nx_full = 512, # number of points in the x direction for full simulation
               Ny_full = 512, # number of points in the y direction for full simulation
               Nz_full = 64,  # number of points in the z direction for full simulation
               nTf = 12,      # simulation length in unit of inertial period
               cfl = 0.9,     # Courant-Friedrichs–Lewy number
               max_Δt = 5minutes, # max time step
               start_from_restratified = false,
               pickup_checkpoint = true,
               ckpdir_affix = "Regular",
               )

    ShortFront = (; commons...,
               Lx = 1kilometers,
               Ly = 250meters,
               Nx_full = 512,
               Ny_full = 128,
               Nz_full = 128,
               nTf = 4,
               cfl = 0.95,
               max_Δt = 0.5minutes,
               start_from_restratified = false,
               pickup_checkpoint = false,
               ckpdir_affix = "ShortFront",
               )

    NoFront = (; commons...,
               Lx = 250meters,
               Ly = 250meters,
               Nx_full = 128,
               Ny_full = 128,
               Nz_full = 128,
               nTf = 4,
               cfl = 0.65,
               max_Δt = 0.5minutes,
               start_from_restratified = false,
               pickup_checkpoint = false,
               ckpdir_affix = "Regular",
               )
end


function decode_casename(casename)
    fres, b_gradient, heat_flux, wind_stress, wind_direction, Stokes_flag = split(casename, "_")
    case_type = fres[1]
    Gb_type   = ifelse(case_type=='n', 'N', 'M')
    coarsen_h = parse(Int64, fres[2])
    coarsen_z = parse(Int64, fres[3])
    Gb = parse(Float64, lstrip(b_gradient, Gb_type)) / 1e8 # [s⁻²] horizontal or vertical buoyancy gradient
    Q₀ = parse(Float64, lstrip(heat_flux,      'Q'))
    τ₀ = parse(Float64, lstrip(wind_stress,    'W')) / 1e3
    if startswith(wind_direction, 'D')
        θ₀ = parse(Float64, lstrip(wind_direction, 'D'))
        wind_oscillation = false
    elseif startswith(wind_direction, 'O')
        θ₀ = parse(Float64, lstrip(wind_direction, 'O'))
        wind_oscillation = true
    else
        throw(DomainError(wind_direction, "wind direction identifier not supoorted"))
    end
    use_Stokes = parse(Bool, Stokes_flag[end])
    return case_type, coarsen_h, coarsen_z, Gb, Q₀, τ₀, θ₀, wind_oscillation, use_Stokes
end


function enrich_parameters(params, casename)
    case_type, coarsen_h, coarsen_z, Gb, Q₀, τ₀, θ₀, wind_oscillation, use_Stokes = decode_casename(casename)
    nTf        = ifelse(θ₀==90, 3, 4)
    Lz         = ifelse((θ₀==90) & (case_type=='s'), 84meters, 147meters) # depth
    sponge_σ   = round(√((Lz/7)^2 / 6 / 2), sigdigits=3) # [m] sponge layer Gaussian mask width (thickness: Lz/7, tapering to e⁻⁶)
    σ_wind     = ifelse(wind_oscillation, params.f, 0)
    wind_waves = ifelse((τ₀==0) & use_Stokes, false, true)

    M²  = ifelse(case_type=='n', 0,   Gb)
    N₀² = ifelse(case_type=='n', Gb,  params.RiB₀*(M²)^2/(params.f)^2) # [s⁻²] mixed layer buoyancy frequency
    N₁² = ifelse(case_type=='n', N₀², params.RiB₁*N₀²) # [s⁻²] thermocline buoyancy frequency
    Tf  = 2π/params.f   # [s] inertial period
    t₀  = ifelse(params.pickup_checkpoint, 6*Tf, 0) # [s] when to introduce wind stress
    tᵣ  = √2*Tf # [s] length of the linear ramp of wind forcing
    ckp_interval = 1*Tf # [s] how often to checkpoint
    ∂v∂z_cgeo    = ifelse(params.counter_geo_stress, M²/params.f, 0) # vertical gradient of along-front velocity at top/bottom

    Nx = params.Nx_full ÷ coarsen_h
    Ny = params.Ny_full ÷ coarsen_h
    Nz = params.Nz_full ÷ coarsen_z

    τ₀ˣ   = τ₀ * cosd(θ₀)
    τ₀ʸ   = τ₀ * sind(θ₀)
    Cd    = 1.2e-3 # neutral drag coefficient, LP81, #Fig. 9 from Edson et al. 2013
    U₁₀   = ifelse(wind_waves, √(τ₀ / params.ρₐ / Cd), 5) # [m s⁻¹] surface wind speed used to specify surface waves
    Uˢ    = ifelse(use_Stokes, 0.068, 0 )#0.0155 * U₁₀, 0) # [m s⁻¹] surface Stokes drift velocity
    #Dˢ    = ifelse(use_Stokes, 0.14 * U₁₀^2 / g_Earth, 1) # [m] e-folding scale of the Stokes drift (avoid blowing up)
    Dˢ    = ifelse(use_Stokes, 4.77, 1)
    B₀    = g_Earth * params.αᵀ * Q₀ / (params.ρ₀ * params.cₚ) # [m² s⁻³] surface buoyancy flux
    ustar = √(τ₀ / params.ρ₀) # [m s⁻¹] friction velocity
    wstar = ∛(B₀ * params.hᵢ) # [m s⁻¹] convective velocity
    Vg    = M² / params.f * max(Lz - params.hᵢ, params.hᵢ) # [m s⁻¹] geostrophic velocity scale

    extra_params = Base.@locals()
    delete!(extra_params, :params)
    delete!(extra_params, :Gb)
    delete!(extra_params, :casename)
    return merge(params, NamedTuple(extra_params))
end
