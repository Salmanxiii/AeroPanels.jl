"""
$(SIGNATURES)

Create and return a tuple of cache arrays for the steady solver.
"""
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
    Fa = zeros(TArray, nSSeg + nCSeg) # Aerodynamic forces
    Ra = zeros(TArray, nSSeg + nCSeg) # Location of aerodynamic forces
    Ra .= model.segmentProps.mid

    return (rVertex=rVertex, vVertex=vVertex, b=b, Γp=Γp, Γw=Γw, Γs=Γs, Fa=Fa, Ra=Ra)
end

"""
$(SIGNATURES)

Calculate the deformed mesh coordinates and the kinematic velocity at the mesh vertices.
"""
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
        vVertex[i] = vVertex[i] - vb_g + vs[i] - cross(ωb_g, rVertex[i] - CG)
    end
    return nothing
end

"""
$(SIGNATURES)

Calculate the normal wash on each panel.
"""
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
    # loop over surfaces
    for (s, nc, ns) in sizes
        # Spanwise Segments
        @batch for i in 1:nc
            for j in 1:ns
                i1, i2 = VertexIndex(s, i, j, sizes), VertexIndex(s, i+1, j, sizes)
                i3, i4 = VertexIndex(s, i, j+1, sizes), VertexIndex(s, i+1, j+1, sizes)
                v1, v2, v3, v4 = vVertex[i1], vVertex[i2], vVertex[i3], vVertex[i4]
                # Velocity at mid-segment (1/4 chord of the spanwise segment)
                v = 0.375*(v1 + v2) + 0.125*(v3 + v4)
                m = SpanSegmentIndex(s, i, j, sizes)
                Fa[m] = ρ * Γs[m] * cross(Fa[m] + v, model.segmentProps.r[m])
            end
        end
        # Trailing Edge segments (force is zero because Γs=0)
        for j in 1:ns 
            m = SpanSegmentIndex(s, nc+1, j, sizes)
            Fa[m] = zero(eltype(Fa))
        end
        # Chordwise Segments
        @batch for i in 1:nc
            for j in 1:ns+1
                i1, i2 = VertexIndex(s, i, j, sizes), VertexIndex(s, i+1, j, sizes)
                v1, v2 = vVertex[i1], vVertex[i2]
                v = 0.25*v1 + 0.75*v2
                m = ChordSegmentIndex(s, i, j, sizes)
                Fa[m+nss] = ρ * Γs[m+nss] * cross(Fa[m+nss] + v, model.segmentProps.r[m+nss])
            end
        end
    end
    return nothing
end

"""
    AeroSolve(vb, model; kwargs...)

Solve the steady aerodynamic problem for a given body velocity `vb`.
Returns a cache containing all circulations and forces.
"""
function AeroSolve(vb::AbstractArray, model::AeroModel;
        ωb=SA[0.0, 0.0, 0.0],
        δc=zeros(eltype(vb), length(model.controlSurfaces)),
        dδc=zeros(eltype(vb), length(model.controlSurfaces)),
        rs=fill(zero(Point3{eltype(vb)}), model.sizes.totalVertices),
        vs=fill(zero(Point3{eltype(vb)}), model.sizes.totalVertices),
        ρ=1.225)

    T = promote_type(eltype(vb), eltype(ωb), eltype(δc), eltype(dδc))
    cache = CreateCacheArrays(model, T)
    AeroSolve!(cache, vb, ωb, δc, dδc, rs, vs, model, ρ)
    return cache
end

function AeroSolve(vb::AbstractArray{A}, ωb::AbstractArray{B},
    δc::AbstractArray{C}, dδc::AbstractArray{D},
    rs::AbstractArray, vs::AbstractArray, model, ρ=1.225) where {A, B, C, D}

    AeroSolve(vb, model; ωb=ωb, δc=δc, dδc=dδc, rs=rs, vs=vs, ρ=ρ)
end

"""
$(SIGNATURES)

In-place version of `AeroSolve`.
"""
function AeroSolve!(cache, vb, ωb, δc, dδc, rs, vs, model, ρ=1.225)
    AddSteadyKinematics!(cache.rVertex, cache.vVertex, vb, ωb, δc, dδc, rs, vs, model)
    CalculateNormalwash!(cache.b, cache.rVertex, cache.vVertex, model)
    SolveCirculation!(cache.Γp, cache.Γw, cache.Γs, cache.b, model)
    CalculateAerodynamicForce!(cache.Fa, cache.Γp, cache.Γw, cache.Γs, cache.vVertex, model, ρ)
    return nothing
end

"""
$(SIGNATURES)

Return the total force and moment at the reference point (CG).
"""
function GetTotalForces(cache, model)
    CG = model.modelProperties.CG
    return IntegrateLoad(cache.Fa, cache.Ra, CG)
end

"""
$(SIGNATURES)

Return the stability axis force and moment coefficients.
"""
function GetStabilityCoefficients(cache, vb, ρ, model)
    mProps = model.modelProperties
    F, M = GetTotalForces(cache, model)
    α, _, _ = AerodynamicAngles(vb)    
    Fstab = GeometryToStabilityAxis(F, α, mProps)
    Mstab = GeometryToStabilityAxis(M, α, mProps)
    CFstab, CMstab = ConvertToCoefficients(Fstab, Mstab, vb, ρ, mProps)
    return CFstab, CMstab
end

"""
$(SIGNATURES)

Calculate aerodynamic loads at monitor points.
"""
function MonitorPointLoads!(Fmp, cache, model::AeroModel)
    for (i, mp) in enumerate(model.monitorPoints)
        indices = mp.segmentIndices
        Fvec = @view cache.Fa[indices]
        Rvec = @view cache.Ra[indices]
        F, M = IntegrateLoad(Fvec, Rvec, mp.origin)
        Fmp[6(i-1)+1:6(i-1)+6] .= [mp.orientation * F; mp.orientation * M]
    end
    return nothing
end
