using Oceananigans
using Oceananigans.Operators
using Oceananigans.Advection: div_Uc
using Oceananigans.TurbulenceClosures: ∇_dot_qᶜ, diffusivity, viscosity 
using Oceananigans.Utils: SumOfArrays

@inline ψ′²(i, j, k, grid, ψ, ψ̄) = @inbounds (ψ[i, j, k] - ψ̄[i, j, k])^2
@inline kinetic_energy_ccc(i, j, k, grid, u, v, w, U, V, W) = (ℑxᶜᵃᵃ(i, j, k, grid, ψ′², u, U) +
                                                               ℑyᵃᶜᵃ(i, j, k, grid, ψ′², v, V) +
                                                               ℑzᵃᵃᶜ(i, j, k, grid, ψ′², w, W)) / 2

function get_output_tuple(model)
    grid = model.grid
    advection = model.advection

    u, v, w = model.velocities
    b = model.tracers.b
    ζ = Field(∂x(v) - ∂y(u))

    # Along-front average
    B  = Field(Average(b, dims=2))
    U  = Field(Average(u, dims=2))
    V  = Field(Average(v, dims=2))
    W  = Field(Average(w, dims=2))
    wb = @at (Center, Center, Center) Field(w * b)
    wbym = Field(Average(wb, dims=2))

    # Bulk N² budget (see Thomas & Ferrari 2008, Eq. 3)
    # note the background advection of background buoyancy is exactly zero
    hrzt_velocities = (u = SumOfArrays{2}(u, model.background_fields.velocities.u),
                       v = SumOfArrays{2}(v, model.background_fields.velocities.v),
                       w = ZFaceField(grid))
    vert_velocities = (u = XFaceField(grid),
                       v = YFaceField(grid),
                       w = SumOfArrays{2}(w, model.background_fields.velocities.w))
    bt = Field(b + model.background_fields.tracers.b)
    badH = KernelFunctionOperation{Center, Center, Center}(div_Uc, grid, advection, hrzt_velocities, bt)
    badV = KernelFunctionOperation{Center, Center, Center}(div_Uc, grid, advection, vert_velocities, bt)
    bdia = KernelFunctionOperation{Center, Center, Center}(∇_dot_qᶜ, grid, model.closure, model.diffusivity_fields, Val(1), b,
                                                           model.clock, fields(model), model.buoyancy)
    BadH = Field(Average(badH, dims=(1,2)))
    BadV = Field(Average(badV, dims=(1,2)))
    Bdia = Field(Average(bdia, dims=(1,2)))
    
    # Equivalence of volume averaged PV through divergence theorem
    # If use Integral, we have to wrap the vorticity and buoyancy into fields first
    vt  = Field(v + model.background_fields.velocities.v)
    btf = @at (Face, Face, Face) bt
    ωaᶻ = @at (Face, Face, Face) Field(f + ζ)
    ωaˣ = @at (Face, Face, Face) Field(∂y(w) - ∂z(vt))
    PVfz = Field(Average(ωaᶻ*btf, dims=(1,2)))
    PVfx = Field(Average(ωaˣ*btf, dims=(2,3)))
    #pv = ErtelPotentialVorticity(model; location=(Face, Face, Face), add_background=true)
    #PV = Field(Average(pv, dims=2))
    #wpv_op = @at (Face, Face, Face) w * pv
    #Jz_adv = Average(wpv_op, dims=2)
    #Jz_fric =
    #Jz_dia
    
    # Horizontal average of stratification 
    b_hrzt_mean = Field(Average(B, dims=1))
    ∂b∂z_bcs = FieldBoundaryConditions(grid, (Nothing, Nothing, Face);
                                       top    = OpenBoundaryCondition(-B₀/κ₀),
                                       bottom = OpenBoundaryCondition(N₁²))
    ∂b∂z_hrzt_mean = Field(∂z(b_hrzt_mean), boundary_conditions=∂b∂z_bcs)
    N²hm = Field(@at (Nothing, Nothing, Center) ∂b∂z_hrzt_mean)
   
    # Turbulent Kinetic Energy
    u_hrzt_mean = Field(Average(U, dims=1))
    v_hrzt_mean = Field(Average(V, dims=1))
    w_hrzt_mean = Field(Average(W, dims=1))
    tke_op = KernelFunctionOperation{Center, Center, Center}(kinetic_energy_ccc,
                                                             grid, u, v, w, u_hrzt_mean, v_hrzt_mean, w_hrzt_mean)
    TKE = Field(Average(tke_op, dims=(1,2,3)))

    # Diffusivity & viscosity
    νₑc = sum(viscosity(model.closure, model.diffusivity_fields))
    κₑc = sum(diffusivity(model.closure, model.diffusivity_fields, Val(:b)))
    νₑ  = Field(@at (Center, Center, Face) νₑc)
    κₑ  = Field(@at (Center, Center, Face) κₑc)

    # Assemble outputs
    fields_slice = Dict("u" => u, "v" => v, "w" => w, "b" => b, "ζ" => ζ, "νₑ" => νₑ, "κₑ" => κₑ)
    line_mean = Dict("B" => B, "U" => U, "V" => V, "W" => W, "wbym" => wbym)
    slice_mean = Dict("N²hm" => N²hm, "PVfz" => PVfz, "PVfx" => PVfx, "BadH" => BadH, "BadV" => BadV, "Bdia" => Bdia)#,
    volume_mean = Dict("TKE" => TKE)#,
    fields_mean = merge(line_mean, slice_mean, volume_mean)
    return fields_slice, fields_mean 
end

#adv_cfl(model) = AdvectiveCFL(simulation.Δt)(model)
#dif_cfl(model) = DiffusiveCFL(simulation.Δt)(model)

#w′ = w - w_hrzt_mean
#b′ = b - b_hrzt_mean
#w′b′_op = @at (Center, Center, Center) w′ * b′
#w′b′ = Average(w′b′_op, dims=(2νₑ = @at (Center, Nothing, Face) model.diffusivity_fields[2].νₑ

#function write_to_ds(dsname, varname, data; mode="c", coords=("xC", "yC", "zC"), dtype=Float64)
#    NCD.NCDataset(dsname, mode) do ds
#        NCD.defVar(ds, varname, data, coords)
#    end
#end

#function save_bak(dsname)
#    Vbak = Field(model.background_fields.velocities.v + v*0)
#    Bbak = Field(model.background_fields.tracers.b + b*0)
#    compute!(Vbak)
#    compute!(Bbak)
#    write_to_ds(dsname, "V", interior(Vbak), mode="c", coords=("xC", "yF", "zC"))
#    write_to_ds(dsname, "B", interior(Bbak), mode="a", coords=("xC", "yC", "zC"))
#end
