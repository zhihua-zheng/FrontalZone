using Oceananigans.Fields: ZeroField
#using Oceananigans.TurbulenceClosures: diffusivity, viscosity 
#using Oceananigans.Models.NonhydrostaticModels: tracer_tendency
#using Oceananigans.Utils: SumOfArrays

#using Oceanostics.FlowDiagnostics: ErtelPotentialVorticity#, RichardsonNumber, RossbyNumber
using Oceanostics.TKEBudgetTerms: PressureRedistributionTerm
#using Oceanostics.PotentialEnergyEquationTerms: PotentialEnergy


fff_scratch = Field{Face, Face, Face}(grid)
ccc_scratch = Field{Center, Center, Center}(grid)
#cnc_scratch = Field{Center, Nothing, Center}(grid)

#function XYSlice(grid, z)
#    i = Colon()
#    j = Colon()
#    ks = searchsortedfirst.(Ref(znodes(grid, Center(), Center(), Face())), z)
#    return ((i, j, k) for k=ks)
#end

#@inline fψ²_plus_gφ²(i, j, k, grid, f, ψ, g, φ) = f(i, j, k, grid, ψ)^2 + g(i, j, k, grid, φ)^2

#function pseudo_isotropic_viscous_dissipation_rate_ccc(i, j, k, grid, u, v, w, p)
#    dudx² = ∂xᶜᶜᶜ(i, j, k, grid, u)^2
#    dvdy² = ∂yᶜᶜᶜ(i, j, k, grid, v)^2
#    dwdz² = ∂zᶜᶜᶜ(i, j, k, grid, w)^2

#    dudy² = ℑxyᶜᶜᵃ(i, j, k, grid, ∂yᶠᶠᶜ, u)^2
#    dvdx² = ℑxyᶜᶜᵃ(i, j, k, grid, ∂xᶠᶠᶜ, v)^2

#    dudz² = ℑxzᶜᵃᶜ(i, j, k, grid, ∂zᶠᶜᶠ, u)^2
#    dwdx² = ℑxzᶜᵃᶜ(i, j, k, grid, ∂xᶠᶜᶠ, w)^2

#    dvdz² = ℑyzᵃᶜᶜ(i, j, k, grid, ∂zᶜᶠᶠ, v)^2
#    dwdy² = ℑyzᵃᶜᶜ(i, j, k, grid, ∂yᶜᶠᶠ, w)^2

    #dudy²_plus_dvdx² = ℑxyᶜᶜᵃ(i, j, k, grid, fψ²_plus_gφ², ∂yᶠᶠᶜ, u, ∂xᶠᶠᶜ, v)
    #dudz²_plus_dwdx² = ℑxzᶜᵃᶜ(i, j, k, grid, fψ²_plus_gφ², ∂zᶠᶜᶠ, u, ∂xᶠᶜᶠ, w)
    #dvdz²_plus_dwdy² = ℑyzᵃᶜᶜ(i, j, k, grid, fψ²_plus_gφ², ∂zᶜᶠᶠ, v, ∂yᶜᶠᶠ, w)

#    ν = _νᶜᶜᶜ(i, j, k, grid, p.closure, p.diffusivity_fields, p.clock)
#    return ν * (dudx² + dvdy² + dwdz² + dudy² + dvdx² + dudz² + dwdx² + dvdz² + dwdy²)
    #return ν * (dudx² + dvdy² + dwdz² + dudy²_plus_dvdx² + dudz²_plus_dwdx² + dvdz²_plus_dwdy²)
#end

#function PseudoIsotropicKineticEnergyDissipationRate(model; U=0, V=0, W=0,
#                                         location = (Center, Center, Center))

#    validate_location(location, "PseudoIsotropicKineticEnergyDissipationRate")
#    validate_dissipative_closure(model.closure)

#    u, v, w = model.velocities

#    parameters = (closure = model.closure,
#                  diffusivity_fields = model.diffusivity_fields,
#                  clock = model.clock)

#    return KernelFunctionOperation{Center, Center, Center}(pseudo_isotropic_viscous_dissipation_rate_ccc, model.grid,
#                                                           (u - U), (v - V), (w - W), parameters)
#end

#@inline function uᵢ∂ᵢtkeᶜᶜᶜ(i, j, k, grid, velocities, tke)
#    u∂x_tke = ℑxᶜᵃᵃ(i, j, k, grid, ψf, velocities.u, ∂xᶠᶜᶜ, tke)
#    v∂y_tke = ℑyᵃᶜᵃ(i, j, k, grid, ψf, velocities.v, ∂yᶜᶠᶜ, tke)
#    w∂z_tke = ℑzᵃᵃᶜ(i, j, k, grid, ψf, velocities.w, ∂zᶜᶜᶠ, tke)
#    return u∂x_tke + v∂y_tke + w∂z_tke
#end
#
#function TurbulenceRedistributionTerm(model::NonhydrostaticModel, tke; velocities=model.velocities)
#    return KernelFunctionOperation{Center, Center, Center}(uᵢ∂ᵢtkeᶜᶜᶜ, model.grid, velocities, tke)
#end

#@inline function ∂ⱼ_uᵢτᵢⱼᶜᶜᶜ(i, j, k, grid, diffusivity_fields, fields, p)
#    τᵢⱼ∂ⱼ_uᵢ_ccc = viscous_dissipation_rate_ccc(i, j, k, grid, diffusivity_fields, fields, p)
#    uᵢ∂ⱼ_τᵢⱼ_ccc = uᵢ∂ⱼ_τᵢⱼᶜᶜᶜ(i, j, k, grid, p.closure, diffusivity_fields, p.clock, fields, p.buoyancy)
#    return uᵢ∂ⱼ_τᵢⱼ_ccc - τᵢⱼ∂ⱼ_uᵢ_ccc
#end
#
#@inline function SubgridscaleRedistributionTerm(model::NonhydrostaticModel; U=ZeroField(), V=ZeroField(), W=ZeroField())
#    mean_velocities = (u=U, v=V, w=W)
#    model_fields = perturbation_fields(model; mean_velocities...)
#    parameters = (; model.closure,
#                  model.clock,
#                  model.buoyancy)
#    return KernelFunctionOperation{Center, Center, Center}(∂ⱼ_uᵢτᵢⱼᶜᶜᶜ, model.grid,
#                                                           model.diffusivity_fields, model_fields, parameters)
#end


@inline function get_output_tuple(model, uˢ, vˢ, Vg, p; extra_outputs=false)
    uE = Field(model.velocities.u - uˢ)
    vE = ifelse(p.use_background_Vg, Field(model.velocities.v - vˢ),
                                     Field(model.velocities.v - vˢ - Vg))
    w  = model.velocities.w
    b  = model.tracers.b
    #c  = model.tracers.c

    wᶠᶠᶠ = @at (Face,   Face,   Face)   w
    wᶜᶜᶜ = @at (Center, Center, Center) w
    uᶜᶜᶜ = @at (Center, Center, Center) uE 
    vᶜᶜᶜ = @at (Center, Center, Center) vE 

    #bt = Field(b + model.background_fields.tracers.b)
    #Ri = Field(RichardsonNumber(model, uE, vE, w, model.tracers.b))
    #Ro = Field(RossbyNumber(model, uE, vE, w, model.coriolis))
    q = Field(ErtelPotentialVorticityFrontalZone(model, uE, vE, w, b, model.coriolis, M²=p.M²), data=fff_scratch.data)

    # Along-front average
    #bym = Field(Average(b, dims=2), data=cnc_scratch.data)
    #cym = Field(Average(c, dims=2), data=cnc_scratch.data)
    #uym = Field(Average(u_ccc, dims=2), data=cnc_scratch.data)
    #vym = Field(Average(v_ccc, dims=2), data=cnc_scratch.data)
    #wym = Field(Average(w_ccc, dims=2), data=cnc_scratch.data)

    #Riym = Field(Average(Ri, dims=2))
    #Roym = Field(Average(Ro, dims=2))
    #PVym = Field(Average(PV, dims=2))
    
    # Diffusivity & viscosity
    #ufrc = -KernelFunctionOperation{Center, Center, Center}(∂ⱼ_τ₁ⱼ,   grid, model.closure, model.diffusivity_fields,
    #                                                        model.clock, fields(model), model.buoyancy)
    #vfrc = -KernelFunctionOperation{Center, Center, Center}(∂ⱼ_τ₂ⱼ,   grid, model.closure, model.diffusivity_fields,
    #                                                        model.clock, fields(model), model.buoyancy)
    #wfrc = -KernelFunctionOperation{Center, Center, Center}(∂ⱼ_τ₃ⱼ,   grid, model.closure, model.diffusivity_fields,
    #                                                        model.clock, fields(model), model.buoyancy)
    #bdia = -KernelFunctionOperation{Center, Center, Center}(∇_dot_qᶜ, grid, model.closure, model.diffusivity_fields, Val(1), b,
    #                                                        model.clock, fields(model), model.buoyancy)
    #νₑ = sum(viscosity(model.closure, model.diffusivity_fields))
    #κₑ = sum(diffusivity(model.closure, model.diffusivity_fields, Val(:b)))
    #νₑ_ccc = viscosity(model.closure, model.diffusivity_fields)
    #κₑ_ccc = diffusivity(model.closure, model.diffusivity_fields, Val(:b))
    #νₑ_ccf = @at (Center, Center, Face) νₑ_ccc
    #κₑ_ccf = @at (Center, Center, Face) κₑ_ccc
    #νₑhm = Field(Average(νₑ_ccf, dims=(1,2)))
    #κₑhm = Field(Average(κₑ_ccf, dims=(1,2)))

    uusgs = Field(Average(XSubgridscaleNormalStress(model), dims=(1,2)))
    vvsgs = Field(Average(YSubgridscaleNormalStress(model), dims=(1,2)))
    wwsgs = Field(Average(ZSubgridscaleNormalStress(model), dims=(1,2)))
    wusgs = Field(Average(XSubgridscaleVerticalMomentumFlux(model),  dims=(1,2)))
    wvsgs = Field(Average(YSubgridscaleVerticalMomentumFlux(model),  dims=(1,2)))
    wbsgs = Field(Average(SubgridscaleVerticalTracerFlux(model, :b), dims=(1,2)))
    #wcsgs = Field(Average(SubgridscaleVerticalTracerFlux(model, :c), dims=(1,2)))

    # Correlations
    uh = Field(Average(uE, dims=(1,2)))
    vh = Field(Average(vE, dims=(1,2)))
    u′ = Field(uE - uh)
    v′ = Field(vE - vh)
    mean_velocities = (u=uh, v=vh, w=ZeroField())
    pert_velocities = (u=u′, v=v′, w=w)
    geo_velocities  = (u=ZeroField(), v=Vg, w=ZeroField())

    uhᶜ = Field(Average(uᶜᶜᶜ, dims=(1,2)))
    vhᶜ = Field(Average(vᶜᶜᶜ, dims=(1,2)))
    bh  = Field(Average(b,    dims=(1,2)))
    #ch  = Field(Average(c,    dims=(1,2)))
    qh  = Field(Average(q,    dims=(1,2)))
    u′ᶜ = Field(uᶜᶜᶜ - uhᶜ)
    v′ᶜ = Field(vᶜᶜᶜ - vhᶜ)
    b′  = Field(b    - bh )
    #c′  = Field(c    - ch )
    q′  = Field(q    - qh )
    uut = Field(Average(u′ᶜ  * u′ᶜ , dims=(1,2)))
    vvt = Field(Average(v′ᶜ  * v′ᶜ , dims=(1,2)))
    wwt = Field(Average(wᶜᶜᶜ * wᶜᶜᶜ, dims=(1,2)))
    bbt = Field(Average(b′   * b′  , dims=(1,2)))
    #cct = Field(Average(c′   * c′  , dims=(1,2)))
    wut = Field(Average(wᶜᶜᶜ * u′ᶜ , dims=(1,2)))
    wvt = Field(Average(wᶜᶜᶜ * v′ᶜ , dims=(1,2)))
    wbt = Field(Average(wᶜᶜᶜ * b′  , dims=(1,2)))
    #wct = Field(Average(wᶜᶜᶜ * c′  , dims=(1,2)))
    wqt = Field(Average(wᶠᶠᶠ * q′  , dims=(1,2)))
    uvt = Field(Average(u′ᶜ  * v′ᶜ , dims=(1,2)))

    #Tsgshm = Field(Average(SubgridscaleRedistributionTerm(model, U=uh, V=vh), dims=(1,2)))
    eps     = Field(KineticEnergyDissipation(model, energy_vel=pert_velocities), data=ccc_scratch.data)
    TKE_eps = Field(Average(eps, dims=(1,2)))
    TKE_prs = Field(Average(PressureRedistributionTerm(model), dims=(1,2)))
    TKE_tur = Field(Average(KineticEnergyAdvection(model, velocities=pert_velocities, energy_vel=pert_velocities), dims=(1,2)))
    TKE_sgs = Field(Average(KineticEnergyStress(model, energy_vel=pert_velocities), dims=(1,2)))
    TKE_spg = Field(Average(KineticEnergyForcing(model, energy_vel=pert_velocities), dims=(1,2)))

    MKE_eps = Field(Average(KineticEnergyDissipation(model, energy_vel=mean_velocities), dims=(1,2)))
    MKE_tur = Field(Average(KineticEnergyAdvection(model, velocities=pert_velocities, energy_vel=mean_velocities), dims=(1,2)))
    MKE_sgs = Field(Average(KineticEnergyStress(model, energy_vel=mean_velocities), dims=(1,2)))
    MKE_spg = Field(Average(KineticEnergyForcing(model, energy_vel=mean_velocities), dims=(1,2)))

    CKE_eps = Field(Average(KineticEnergyDissipation(model, energy_vel=geo_velocities), dims=(1,2)))
    CKE_tur = Field(Average(KineticEnergyAdvection(model, velocities=pert_velocities, energy_vel=geo_velocities), dims=(1,2)))
    CKE_sgs = Field(Average(KineticEnergyStress(model, energy_vel=geo_velocities), dims=(1,2)))
    CKE_spg = Field(Average(KineticEnergyForcing(model, energy_vel=geo_velocities), dims=(1,2)))

    # Surface fluxes
    Qu = Field(Average(SurfaceMomentumFlux(model, :u), dims=(1,2)))
    Qv = Field(Average(SurfaceMomentumFlux(model, :v), dims=(1,2)))
    Qb = Field(Average(SurfaceTracerFlux(model,   :b), dims=(1,2)))

    # Equivalence of volume averaged PV through divergence theorem
    # If use Integral, we have to wrap the vorticity and buoyancy into fields first
    #Bbkf = Field{Face, Face, Face}(grid)
    #set!(Bbkf, (x, y, z) -> (-pm.M² * x))
#    vfcf = @at (Face, Center, Face) v
#    wffc = @at (Face, Face, Center) w
#    ωᶻ   = @at (Center, Center, Center) Field(∂x(v) - ∂y(u))
#    ωˣ   = Field(∂y(wffc) - ∂z(vfcf))
#    PVfz = Field(Average(ωᶻ*bt, dims=(1,2)))
#    RVx  = Field(Average(ωˣ,    dims=(2,3)))
    #pv = ErtelPotentialVorticity(model; location=(Face, Face, Face), add_background=true)
    #PV = Field(Average(pv, dims=2))
    #wpv_op = @at (Face, Face, Face) w * pv
    #Jz_adv = Average(wpv_op, dims=2)
    #Jz_fric =
    #Jz_dia

    # Assemble outputs
    fields_slice = Dict("u" => uᶜᶜᶜ, "v" => vᶜᶜᶜ, "w" => wᶜᶜᶜ, "b" => b, "q" => q, "eps" => eps)#, "c" => c)
#                        "ufrc" => ufrc, "vfrc" => vfrc, "wfrc" => wfrc, "bdia" => bdia, "νₑ" => νₑ, "κₑ" => κₑ)
    fields_mean = Dict(
#                       "bym" => bym, "cym" => cym, "uym" => uym, "vym" => vym, "wym" => wym,
                       "u" => uhᶜ, "v" => vhᶜ, "b" => bh, "q" => qh,# "c" => ch,
                       "uut" => uut, "vvt" => vvt, "wwt" => wwt, "bbt" => bbt,# "cct" => cct,
                       "wut" => wut, "wvt" => wvt, "uvt" => uvt, "wbt" => wbt, "wqt" => wqt,# "wct" => wct,
                       "uusgs" => uusgs, "vvsgs" => vvsgs, "wwsgs" => wwsgs,
                       "wusgs" => wusgs, "wvsgs" => wvsgs, "wbsgs" => wbsgs,# "wcsgs" => wcsgs,
                       "Qu" => Qu, "Qv" => Qv, "Qb" => Qb,
                       "TKE_eps" => TKE_eps, "TKE_tur" => TKE_tur, "TKE_sgs" => TKE_sgs, "TKE_spg" => TKE_spg, "TKE_prs" => TKE_prs,
                       "MKE_eps" => MKE_eps, "MKE_tur" => MKE_tur, "MKE_sgs" => MKE_sgs, "MKE_spg" => MKE_spg,
                       "CKE_eps" => CKE_eps, "CKE_tur" => CKE_tur, "CKE_sgs" => CKE_sgs, "CKE_spg" => CKE_spg)
#                       "Roym" => Roym, "Riym" => Riym, "PVym" => PVym,
#                       "PVfz" => PVfz, "RVx" => RVx)

    if extra_outputs
        pHSA = Field(Average(model.pressures.pHY′, dims=2))
        pNHS = Field(Average(model.pressures.pNHS, dims=2))

        # Bulk N² budget (see Thomas & Ferrari 2008, Eq. 3)
        # note the background advection of background buoyancy is exactly zero
        total_velocities = (u = SumOfArrays{2}(u, model.background_fields.velocities.u),
                            v = SumOfArrays{2}(v, model.background_fields.velocities.v),
                            w = SumOfArrays{2}(w, model.background_fields.velocities.w))
        badv = -KernelFunctionOperation{Center, Center, Center}(div_Uc, grid, advection, total_velocities, b)
        badb = -KernelFunctionOperation{Center, Center, Center}(div_Uc, grid, advection, velocities, model.background_fields.tracers.b)
        bdia = -KernelFunctionOperation{Center, Center, Center}(∇_dot_qᶜ, grid, model.closure, model.diffusivity_fields, Val(1), b,
                                                                model.clock, fields(model), model.buoyancy)
        ∂ₜb = KernelFunctionOperation{Center, Center, Center}(tracer_tendency, grid, Val(1), Val(:b), advection, model.closure,
                                                              b.boundary_conditions.immersed, model.buoyancy, model.biogeochemistry,
                                                              model.background_fields, velocities, model.tracers, model.auxiliary_fields,
                                                              model.diffusivity_fields, model.forcing.b, model.clock)
        Badv = Field(Average(badv, dims=(1,2)))
        Badb = Field(Average(badb, dims=(1,2)))
        Bdia = Field(Average(bdia, dims=(1,2)))
        ∂ₜB  = Field(Average(∂ₜb,  dims=(1,2)))

        fields_mean_extra = Dict("pHSA" => pHSA, "pNHS" => pNHS, "Badv" => Badv, "Badb" => Badb, "Bdia" => Bdia, "∂ₜB" => ∂ₜB)
        fields_mean = merge(fields_mean, fields_mean_extra)
    end
    return fields_slice, fields_mean
end
