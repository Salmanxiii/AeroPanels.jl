# vb: body velocity, ab: acceleration, ωb: angular velocity, dωb: angular acceleration,
# δc: control surface deflection, dδc: CS deflection derivative, ddδc: CS deflection double derivaitve,
# vd: disturbance velocity at mesh vertices, ad: disturbance acceleration

function CreateCacheArrays(model, T)
    nVert = model.sizes.totalVertices
    nPan = model.sizes.totalPanels
    nWake = sum(ns for (s, nc, ns) in model.sizes)
    nSSeg = model.sizes.totalSpanSegments
    nCSeg = model.sizes.totalChordSegments

    TArray = Point3{T}
    rVertex = zeros(TArray, nVert)
    vVertex = zeros(TArray, nVert)
    b = zeros(T, nPan)
    Γp = zeros(T, nPan)
    Γw = zeros(T, nWake)
    Γs = zeros(T, nSSeg + nCSeg)
    Fa = zeros(TArray, nSSeg + nCSeg)

    return (rVertex=rVertex, vVertex=vVertex, b=b, Γp=Γp, Γw=Γw, Γs=Γs, Fa=Fa)
end

function AddSteadyKinematics!(rVertex, vVertex,
        vb, ωb, δc, dδc, rs, vs, model)

    rVertex .= coordinates(model.mesh)
    fill!(vVertex, zero(eltype(vVertex)))
    CSDisplacement!(rVertex, δc, model)
    CSVelocity!(vVertex, dδc, model)

    bodyAxis = model.modelProperties.bodyFixedCS # body axis coordinate
    CG = model.modelProperties.CG
    # velocities in geometry frame
    (vb_g, ωb_g) = (bodyAxis' * component for component in (vb, ωb))

    for i in eachindex(rVertex)
        rVertex[i] = rVertex[i] + rs[i]
        vVertex[i] = vVertex[i] + vb_g + vs[i] + cross(ωb_g, rVertex[i] - CG)
    end
    return nothing
end

function CalculateNormalwash!(b, rVertex, vVertex, model)
    # Normalwash calculation: calculate new normals and interpolate velocity at collocation point
    for (s, nc, ns) in model.sizes
        @batch for i in 1:nc
            for j in 1:ns
                p = PanelIndex(s, i, j, model.sizes)
                i1, i2, i3, i4 = PanelVertexIndices(s, i, j, model.sizes)
                r1, r2, r3, r4 = rVertex[i1], rVertex[i2], rVertex[i3], rVertex[i4]
                v1, v2, v3, v4 = vVertex[i1], vVertex[i2], vVertex[i3], vVertex[i4]

                normalUnitVec = normalize(cross(r3 - r1, r2 - r4))
                v = 0.125*(v1 + v2) + 0.375*(v3 + v4)
                b[p] = - dot(v, normalUnitVec)
            end
        end
    end
    return nothing
end


"""
$(SIGNATURES)

Solve for panel, wake, and segment circulations in-place.
"""
function SolveCirculation!(Γp, Γw, Γs, b, model::AeroModel)
    # Solve for panel circulation
    ldiv!(Γp, model.AIC, b)
    
    # Extract wake circulation
    index = TEPanelIndex(model.sizes)
    Γw .= @view Γp[index]
    
    # Calculate segment circulation
    SegmentCirculation!(Γs, Γp, model.segmentProps)
    return nothing
end

"""
$(SIGNATURES)

Calculate aerodynamic forces on segments in-place into vector `Fa`.
"""
function CalculateAerodynamicForce!(Fa, Γp, Γw, Γs, vVertex, model::AeroModel, ρ)
    # Calculate induced velocity at segments and storing it temporarily in Fa
    # Non-allocating version of equation v = AIC3r * Γr + AIC3w * Γw
    Fa .= model.segmentProps.aic3Ring * Γp
    mul!(Fa, model.segmentProps.aic3Wake, Γw, 1.0, 1.0)
    
    nss = model.segmentProps.nSpanSegments
    sizes = model.sizes
    # loop over Spanwise Segments
    for (s, nc, ns) in model.sizes
        @batch for i in 1:nc
            for j in 1:ns
                i1, i2, i3, i4 = PanelVertexIndices(s, i, j, sizes)
                v1, v2, v3, v4 = vVertex[i1], vVertex[i2], vVertex[i3], vVertex[i4]
                v = 0.375*(v1 + v2) + 0.125*(v3 + v4)
                m = SpanSegmentIndex(s, i, j, sizes)
                Fa[m] = ρ * Γs[m] * cross(Fa[m] + v, model.segmentProps.r[m])
            end
            for k in 1:ns+1
                v1 = vVertex[VertexIndex(s, i, k, sizes)]
                v2 = vVertex[VertexIndex(s, i+1, k, sizes)]
                v = 0.25*v1 + 0.75*v2
                Fa[m+nss] = ρ * Γs[m+nss] * cross(Fa[m+nss] + v, model.segmentProps.r[m+nss])
            end
        end
        # for i = nc+1 trailing edge segments force = 0 because Γs = 0
        for j in 1:ns 
                m = SpanSegmentIndex(s, nc+1, j, sizes)
                Fa[m] = zeros(eltype(Fa))
        end
    end
    return nothing
end

function AeroSolve(vb::AbstractArray{A}, ωb::AbstractArray{B},
        δc::AbstractArray{C}, dδc::AbstractArray{D},
        rs::AbstractArray, vs::AbstractArray, model, ρ=1.225)

    T = promote_type(A, B, C, D)
    cache = CreateCacheArrays(model, T)
    AeroSolve!(cache, vb, ωb, δc, dδc, rs, vs, model, ρ)
    return cache
end

function AeroSolve!(cache, vb, ωb, δc, dδc, rs, vs, model, ρ=1.225)
    AddSteadyKinematics!(cache.rVertex, cache.vVertex, vb, ωb, δc, dδc, rs, vs, model)
    CalculateNormalwash!(cache.b, cache.rVertex, cache.vVertex, model)
    SolveCirculation!(cache.Γp, cache.Γw, cache.Γs, cache.b, model)
    CalculateAerodynamicForce!(cache.Fa, cache.Γp, cache.Γw, cache.Γs, cache.vVertex, model, ρ)
    return nothing
end

function GetTotalForces(Fa, model)
    F = sum(Fa)
    r = AerodynamicLoadLocation(model)
    CG = model.modelProperties.CG
    M = sum(cross(r[i]-CG, Fa[i]) for i in eachindex(r, Fa))
    return F, M
end

function GetCoefficients(Fa, vb, ρ,model)
    mProps = model.modelProperties

    F, M = GetTotalForces(Fa, model)
    α, β, V = AerodynamicAngles(vb)
    Q = 0.5 * ρ * V^2
    
    CFstab = BodyFixedToStabilityAxis(F, α) / Q / mProps.S
    CMstab = (GeometryToStabilityAxis(M, α, mProps) / Q / mProps.S) ./ SA[mProps.b, mProps.c, mProps.b] 
    return CFstab, CMstab
end