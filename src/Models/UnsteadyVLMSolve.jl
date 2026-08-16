# src/Models/UnsteadyVLMSolve.jl

function AddUnsteadyKinematics!(aVertex, rVertex,
    ab, dωb, ddδc, as, model)

    fill!(aVertex, zero(eltype(aVertex)))
    CSAcceleration!(aVertex, ddδc, model)

    rRef = model.modelProperties.rRef
    bodyAxis = model.modelProperties.bodyFixedCS
    (ab_g, dωb_g) = (bodyAxis' * component for component in (ab, dωb))

    for i in eachindex(aVertex)
        aVertex[i] = aVertex[i] - ab_g + as[i] - cross(dωb_g, rVertex[i] - rRef)
    end
    return nothing
end

function CalculateDNormalwash!(db, rVertex, vVertex, aVertex, model)
    for (s, nc, ns) in model.boundMesh.sizes
        @batch for i in 1:nc
            for j in 1:ns
                p = PanelIndex(s, i, j, model.boundMesh.sizes)
                i1, i2, i3, i4 = PanelVertexIndices(s, i, j, model.boundMesh.sizes)
                r1, r2, r3, r4 = rVertex[i1], rVertex[i2], rVertex[i3], rVertex[i4]
                v1, v2, v3, v4 = vVertex[i1], vVertex[i2], vVertex[i3], vVertex[i4]
                a1, a2, a3, a4 = aVertex[i1], aVertex[i2], aVertex[i3], aVertex[i4]

                normalVec = cross(r3 - r1, r2 - r4)
                normalNorm = sqrt(dot(normalVec, normalVec))
                normalUnitVec = normalVec / normalNorm
                
                v = 0.125*(v1 + v2) + 0.375*(v3 + v4)
                a = 0.125*(a1 + a2) + 0.375*(a3 + a4)

                dNormalVec = cross(v3 - v1, r2 - r4) + cross(r3 - r1, v2 - v4)
                dNormalUnitVec = (dNormalVec - normalUnitVec * dot(normalUnitVec, dNormalVec)) / normalNorm
                db[p] = -(dot(a, normalUnitVec) + dot(v, dNormalUnitVec))
            end
        end
    end
    return nothing
end

function SolveUnsteadyCirculation!(dΓb, Γb, Γw, Γs, b, db, V, model::UnsteadyAeroModel2D, Γw1)
    mul!(Γb, model.L3, Γw1)
    mul!(Γb, model.L4, b, -1.0, 1.0)

    mul!(dΓb, model.L5, Γw1)
    mul!(dΓb, model.L6, b, 1.0, 1.0)
    dΓb .*= V
    mul!(dΓb, model.L4, db, -1.0, 1.0)

    Γw[model.Γw1Indices] .= Γw1
    Γw[model.Γw0Indices] .= @view Γb[model.ΓbTEIndices]

    SegmentCirculation!(Γs, Γb, model.segmentProps)
    return nothing
end

function CalculateUnsteadyAerodynamicForce!(Fa, dΓb, Γb, Γw, Γs, vVertex, model::UnsteadyAeroModel2D, ρ)
    nss = model.boundMesh.sizes.totalSpanSegments
    ncs = model.boundMesh.sizes.totalChordSegments
    nts = nss + ncs
    
    FaSteady = @view Fa[1:nts]
    CalculateAerodynamicForce!(FaSteady, Γb, Γw, Γs, vVertex, model, ρ)

    Fa_uns = @view Fa[nts+1:end]
    Fa_uns .= ρ .* dΓb .* GetPanelAreas(model.boundMesh) .* GetNormalVectors(model.boundMesh)
    return nothing
end

function AeroPanels.InitializeFMU!(array::AbstractFMUArray, cache, model::UnsteadyAeroModel2D, t::Float64; start_from_trim::Bool=false)
    T_el = eltype(array.vb)
    rs = reinterpret(Point3{T_el}, array.rs)
    vs = reinterpret(Point3{T_el}, array.vs)
    as = reinterpret(Point3{T_el}, array.as)
    dδc = Vector{T_el}(array.du_cs)
    ddδc = Vector{T_el}(array.ddu_cs)

    vd = hasproperty(array, :vd) ? SVector{3}(array.vd) : SA[zero(T_el), zero(T_el), zero(T_el)]
    wd = hasproperty(array, :wd) ? SVector{3}(array.wd) : SA[zero(T_el), zero(T_el), zero(T_el)]
    
    v_eff = SVector{3}(array.vb) + vd
    w_eff = SVector{3}(array.wb) + wd

    AddSteadyKinematics!(cache.rVertex, cache.vVertex, v_eff, w_eff, array.u_cs, dδc, rs, vs, model)
    if start_from_trim
        array.circ_w1 .= -model.K8 \ (model.K9 * cache.b)
    else
        fill!(array.circ_w1, zero(T_el))
    end
    EvaluateOutputs!(array, cache, model, t)
    return nothing
end

function AeroPanels.InitializeSteadyState!(array::AbstractFMUArray, cache, model::UnsteadyAeroModel2D, t::Float64)
    InitializeFMU!(array, cache, model, t; start_from_trim=true)
    return nothing
end

function AeroPanels.EvaluateDerivatives!(du::AbstractVector, array::AbstractFMUArray, cache::UnsteadyAeroModel2DCache, model::UnsteadyAeroModel2D, t::Float64)
    V = sqrt(dot(array.vb, array.vb))
    fill!(du, zero(eltype(du)))
    mul!(du, model.K8, array.circ_w1)
    mul!(du, model.K9, cache.b, 1.0, 1.0)
    du .*= V
    return nothing
end

function AeroPanels.EvaluateDerivatives!(array::AbstractFMUArray, cache::UnsteadyAeroModel2DCache, model::UnsteadyAeroModel2D, t::Float64)
    V = sqrt(dot(array.vb, array.vb))
    fill!(array.dcirc_w1, zero(eltype(array.dcirc_w1)))
    mul!(array.dcirc_w1, model.K8, array.circ_w1)
    mul!(array.dcirc_w1, model.K9, cache.b, 1.0, 1.0)
    array.dcirc_w1 .*= V
    return nothing
end

@inline function as_points(arr::AbstractVector{T}) where T
    if isbitstype(T)
        return reinterpret(Point3{T}, arr)
    else
        N = length(arr) ÷ 3
        pts = Vector{Point3{T}}(undef, N)
        @inbounds for i in 1:N
            pts[i] = Point3{T}(arr[3i-2], arr[3i-1], arr[3i])
        end
        return pts
    end
end

function AeroPanels.EvaluateOutputs!(array::AbstractFMUArray, cache::UnsteadyAeroModel2DCache, model::UnsteadyAeroModel2D, t::Float64)
    T_el = eltype(array.vb)
    rs = as_points(array.rs)
    vs = as_points(array.vs)
    as = as_points(array.as)
    dδc = Vector{T_el}(array.du_cs)
    ddδc = Vector{T_el}(array.ddu_cs)

    vd = hasproperty(array, :vd) ? SVector{3}(array.vd) : SA[zero(T_el), zero(T_el), zero(T_el)]
    wd = hasproperty(array, :wd) ? SVector{3}(array.wd) : SA[zero(T_el), zero(T_el), zero(T_el)]
    ad = hasproperty(array, :ad) ? SVector{3}(array.ad) : SA[zero(T_el), zero(T_el), zero(T_el)]
    dwd = hasproperty(array, :dwd) ? SVector{3}(array.dwd) : SA[zero(T_el), zero(T_el), zero(T_el)]

    v_eff = SVector{3}(array.vb) + vd
    w_eff = SVector{3}(array.wb) + wd
    a_eff = SVector{3}(array.ab) + ad
    dw_eff = SVector{3}(array.dwb) + dwd

    AddSteadyKinematics!(cache.rVertex, cache.vVertex, v_eff, w_eff, array.u_cs, dδc, rs, vs, model)
    CalculateNormalwash!(cache.b, cache.rVertex, cache.vVertex, model)

    AddUnsteadyKinematics!(cache.aVertex, cache.rVertex, a_eff, dw_eff, ddδc, as, model)
    CalculateDNormalwash!(cache.db, cache.rVertex, cache.vVertex, cache.aVertex, model)
    
    V = sqrt(dot(array.vb, array.vb))
    SolveUnsteadyCirculation!(cache.dΓb, cache.Γb, cache.Γw, cache.Γs, cache.b, cache.db, V, model, array.circ_w1)
    
    CalculateUnsteadyAerodynamicForce!(cache.Fa, cache.dΓb, cache.Γb, cache.Γw, cache.Γs, cache.vVertex, model, array.rho[1])
    
    Fg, Mg = GetTotalForces(cache, model)
    
    nss = model.boundMesh.sizes.totalSpanSegments
    ncs = model.boundMesh.sizes.totalChordSegments
    nts = nss + ncs
    Fg_unsteady = sum(cache.Fa[nts+1:end])
    
    mProps = model.modelProperties
    Fb = GeometryToBodyAxis(Fg, mProps)
    Mb = GeometryToBodyAxis(Mg, mProps)
    Fb_unsteady = GeometryToBodyAxis(Fg_unsteady, mProps)
    
    array.Fb .= Fb
    array.Mb .= Mb
    array.Fb_unsteady .= Fb_unsteady
    
    vb = SVector{3}(array.vb)
    ρ = array.rho[1]
    α, _, _ = AerodynamicAngles(vb)
    
    CFbody, CMbody = ConvertToCoefficients(Fb, Mb, vb, ρ, mProps)
    array.Cb .= (CFbody..., CMbody...)
    
    Fstab = GeometryToStabilityAxis(Fg, α, mProps)
    Mstab = GeometryToStabilityAxis(Mg, α, mProps)
    CFstab, CMstab = ConvertToCoefficients(Fstab, Mstab, vb, ρ, mProps)
    array.Cs .= (CFstab..., CMstab...)
    
    if length(model.monitorPoints) > 0
        MonitorPointLoads!(array.Fmp, cache, model)
    end
    
    for (i, fvec) in enumerate(cache.Fa)
        array.Fa[3*(i-1)+1] = fvec[1]
        array.Fa[3*(i-1)+2] = fvec[2]
        array.Fa[3*(i-1)+3] = fvec[3]
    end
    
    return nothing
end

function MonitorPointLoads!(Fmp, cache, model::UnsteadyAeroModel2D)
    for (i, mp) in enumerate(model.monitorPoints)
        indices = [mp.segmentIndices; mp.panelIndices]
        Fvec = @view cache.Fa[indices]
        Rvec = @view cache.Ra[indices]
        F, M = IntegrateLoad(Fvec, Rvec, mp.origin)
        Fmp[6(i-1)+1:6(i-1)+6] .= [mp.orientation * F; mp.orientation * M]
    end
    return nothing
end

function AeroPanels.ManualDirectionalDerivative(
    model::UnsteadyAeroModel2D{T},
    array::AbstractFMUArray,
    outputRefs::AbstractVector{<:Integer},
    inputRefs::AbstractVector{<:Integer},
    inputSeed::AbstractVector;
    cache=nothing
) where T
    nStates = NumberOfStates(model)
    
    # Case 1: Pure State-to-StateDerivative override (dcirc_w1 = V * K8 * circ_w1)
    if nStates > 0 && all(ref -> 1 <= ref <= nStates, inputRefs) && all(ref -> (nStates + 1) <= ref <= 2 * nStates, outputRefs)
        V_val = Float64(ForwardDiff.value(norm(array.vb)))
        
        # Build sparse input seed vector for states
        seed_vec = zeros(T, nStates)
        for i in 1:length(inputRefs)
            seed_vec[inputRefs[i]] = T(inputSeed[i])
        end
        
        # Fast sparse matrix-vector product: δ(dcirc_w1) = V * (K8 * δcirc_w1)
        d_out = V_val .* (model.K8 * seed_vec)
        
        res = Vector{T}(undef, length(outputRefs))
        for k in 1:length(outputRefs)
            deriv_idx = outputRefs[k] - nStates
            res[k] = d_out[deriv_idx]
        end
        return res
    end
    
    return nothing
end