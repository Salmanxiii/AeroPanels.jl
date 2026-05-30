"""
$(SIGNATURES)

Create and return a tuple of cache arrays for the unsteady solver.
"""
function CreateCacheArrays(model::UnsteadyAeroModel2D, T)
    nVert = model.sizes.totalVertices
    nPan = model.sizes.totalPanels
    nWake1 = NumberOfStates(model)
    nWake = model.wakeSizes.totalPanels
    nSSeg = model.sizes.totalSpanSegments
    nCSeg = model.sizes.totalChordSegments

    TArray = Point3{T}
    rVertex = zeros(TArray, nVert)
    vVertex = zeros(TArray, nVert)
    aVertex = zeros(TArray, nVert)
    b = zeros(T, nPan)
    db = zeros(T, nPan)
    Γw1 = zeros(T, nWake1)
    dΓw1 = zeros(T, nWake1)
    Γb = zeros(T, nPan)
    dΓb = zeros(T, nPan)
    Γw = zeros(T, nWake)
    Γs = zeros(T, nSSeg + nCSeg)
    Fa = zeros(TArray, nSSeg + nCSeg + nPan)
    Ra = zeros(TArray, nSSeg + nCSeg + nPan)
    Ra[1:nSSeg + nCSeg] .= model.segmentProps.mid
    Ra[nSSeg + nCSeg + 1:end] .= model.panelProperties.rMid

    return (rVertex=rVertex, vVertex=vVertex, aVertex=aVertex, b=b, db=db, 
            Γw1=Γw1, dΓw1=dΓw1, Γb=Γb, dΓb=dΓb, Γw=Γw, Γs=Γs, Fa=Fa, Ra=Ra)
end

"""
$(SIGNATURES)

Calculate the unsteady kinematics (vertex accelerations).
"""
function AddUnsteadyKinematics!(aVertex, rVertex,
    ab, dωb, ddδc, as, model)

    fill!(aVertex, zero(eltype(aVertex)))
    CSAcceleration!(aVertex, ddδc, model)

    bodyAxis = model.modelProperties.bodyFixedCS
    CG = model.modelProperties.CG
    (ab_g, dωb_g) = (bodyAxis' * component for component in (ab, dωb))

    for i in eachindex(aVertex)
        aVertex[i] = aVertex[i] - ab_g + as[i] - cross(dωb_g, rVertex[i] - CG)
    end
    return nothing
end

"""
$(SIGNATURES)

Calculate the time derivative of the normal wash vector `db`.
"""
function CalculateDNormalwash!(db, rVertex, vVertex, aVertex, model)
    # db = - (a . n + v . dn)
    for (s, nc, ns) in model.sizes
        @batch for i in 1:nc
            for j in 1:ns
                p = PanelIndex(s, i, j, model.sizes)
                i1, i2, i3, i4 = PanelVertexIndices(s, i, j, model.sizes)
                r1, r2, r3, r4 = rVertex[i1], rVertex[i2], rVertex[i3], rVertex[i4]
                v1, v2, v3, v4 = vVertex[i1], vVertex[i2], vVertex[i3], vVertex[i4]
                a1, a2, a3, a4 = aVertex[i1], aVertex[i2], aVertex[i3], aVertex[i4]

                normalVec = cross(r3 - r1, r2 - r4)
                normalNorm = norm(normalVec)
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

"""
$(SIGNATURES)

Solve for the time derivatives and intermediate circulations of the unsteady system.
"""
function SolveUnsteadyCirculation!(dΓw1, Γw1, dΓb, Γb, Γw, Γs, b, db, model::UnsteadyAeroModel2D)
    # dΓw1 = K8*Γw1 + K9*b
    mul!(dΓw1, model.K8, Γw1)
    mul!(dΓw1, model.K9, b, 1.0, 1.0)

    # Γb = L3*Γw1 - L4*b
    mul!(Γb, model.L3, Γw1)
    mul!(Γb, model.L4, b, -1.0, 1.0)

    # dΓb = L5*Γw1 + L6*b - L4*db
    mul!(dΓb, model.L5, Γw1)
    mul!(dΓb, model.L6, b, 1.0, 1.0)
    mul!(dΓb, model.L4, db, -1.0, 1.0)

    # Γw = L9 * Γw1 + L10 * b
    mul!(Γw, model.L9, Γw1)
    mul!(Γw, model.L10, b, 1.0, 1.0)

    SegmentCirculation!(Γs, Γb, model.segmentProps)
    return nothing
end

"""
$(SIGNATURES)

Calculate aerodynamic forces for the unsteady model.
"""
function CalculateUnsteadyAerodynamicForce!(Fa, dΓb, Γb, Γw, Γs, vVertex, model::UnsteadyAeroModel2D, ρ)
    nss = model.sizes.totalSpanSegments
    ncs = model.sizes.totalChordSegments
    nts = nss + ncs
    n = model.sizes.totalPanels
    
    FaSteady = @view Fa[1:nts]
    CalculateAerodynamicForce!(FaSteady, Γb, Γw, Γs, vVertex, model, ρ)

    # Unsteady contribution: ρ * dΓ/dt * Area * Normal
    @batch for i in 1:n
        Fa[nts+i] = ρ * dΓb[i] * model.panelProperties.area[i] * model.panelProperties.normal[i]
    end
    return nothing
end

"""
$(SIGNATURES)

In-place version of the unsteady `AeroSolve`.
"""
function AeroSolve!(cache, vb, ab, ωb, dωb, δc, dδc, ddδc, rs, vs, as, Γw1, model::UnsteadyAeroModel2D, ρ=1.225)
    AddSteadyKinematics!(cache.rVertex, cache.vVertex, vb, ωb, δc, dδc, rs, vs, model)
    AddUnsteadyKinematics!(cache.aVertex, cache.rVertex, ab, dωb, ddδc, as, model)
    
    CalculateNormalwash!(cache.b, cache.rVertex, cache.vVertex, model)
    CalculateDNormalwash!(cache.db, cache.rVertex, cache.vVertex, cache.aVertex, model)
    
    SolveUnsteadyCirculation!(cache.dΓw1, Γw1, cache.dΓb, cache.Γb, cache.Γw, cache.Γs, cache.b, cache.db, model)
    CalculateUnsteadyAerodynamicForce!(cache.Fa, cache.dΓb, cache.Γb, cache.Γw, cache.Γs, cache.vVertex, model, ρ)
    return nothing
end

"""
    AeroSolve(vb, ab, Γw1, model; kwargs...)

Solve the unsteady aerodynamic problem for given kinematics and states.
"""
function AeroSolve(vb, ab, Γw1, model::UnsteadyAeroModel2D;
    ωb=SA[0.0, 0.0, 0.0],
    dωb=SA[0.0, 0.0, 0.0],
    δc=zeros(eltype(vb), length(model.controlSurfaces)),
    dδc=zeros(eltype(vb), length(model.controlSurfaces)),
    ddδc=zeros(eltype(vb), length(model.controlSurfaces)),
    rs=fill(zero(Point3{eltype(vb)}), model.sizes.totalVertices),
    vs=fill(zero(Point3{eltype(vb)}), model.sizes.totalVertices),
    as=fill(zero(Point3{eltype(vb)}), model.sizes.totalVertices),
    ρ=1.225)

    T = promote_type(eltype(vb), eltype(ab), eltype(ωb), eltype(dωb), eltype(δc), eltype(dδc), eltype(ddδc), eltype(Γw1))
    cache = CreateCacheArrays(model, T)
    AeroSolve!(cache, vb, ab, ωb, dωb, δc, dδc, ddδc, rs, vs, as, Γw1, model, ρ)
    return cache
end
