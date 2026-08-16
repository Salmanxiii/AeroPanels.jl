# src/Models/SteadyVLMSolve.jl

"""
$(SIGNATURES)

Create and return a tuple of cache arrays for the steady solver.
"""
function CreateCacheArrays(model, T)
    nVert = model.boundMesh.sizes.totalVertices
    nPan = model.boundMesh.sizes.totalPanels
    nWake = model.wakeMesh.wakeSizes.totalPanels
    nSSeg = model.boundMesh.sizes.totalSpanSegments
    nCSeg = model.boundMesh.sizes.totalChordSegments

    TArray = Point3{T}
    rVertex = zeros(TArray, nVert)
    vVertex = zeros(TArray, nVert)
    b = zeros(T, nPan)
    Γp = zeros(T, nPan)
    Γw = zeros(T, nWake)
    Γs = zeros(T, nSSeg + nCSeg)
    Fa = zeros(TArray, nSSeg + nCSeg)
    Ra = zeros(TArray, nSSeg + nCSeg)
    Ra .= model.segmentProps.mid

    return (rVertex=rVertex, vVertex=vVertex, b=b, Γp=Γp, Γw=Γw, Γs=Γs, Fa=Fa, Ra=Ra)
end

"""
$(SIGNATURES)

Calculate the deformed mesh coordinates and the kinematic velocity at the mesh vertices.
"""
function AddSteadyKinematics!(rVertex, vVertex,
        vb, ωb, δc, dδc, rs, vs, model)

    rVertex .= coordinates(GetBoundMesh(model))
    fill!(vVertex, zero(eltype(vVertex)))
    CSDisplacement!(rVertex, δc, model)
    CSVelocity!(vVertex, dδc, model)

    bodyAxis = model.modelProperties.bodyFixedCS
    rRef = model.modelProperties.rRef
    (vb_g, ωb_g) = (bodyAxis' * component for component in (vb, ωb))

    for i in eachindex(rVertex)
        rVertex[i] = rVertex[i] + rs[i]
        vVertex[i] = vVertex[i] - vb_g + vs[i] - cross(ωb_g, rVertex[i] - rRef)
    end
    return nothing
end

"""
$(SIGNATURES)

Calculate the normal wash on each panel.
"""
function CalculateNormalwash!(b, rVertex, vVertex, model)
    for (s, nc, ns) in model.boundMesh.sizes
        @batch for i in 1:nc
            for j in 1:ns
                p = PanelIndex(s, i, j, model.boundMesh.sizes)
                i1, i2, i3, i4 = PanelVertexIndices(s, i, j, model.boundMesh.sizes)
                r1, r2, r3, r4 = rVertex[i1], rVertex[i2], rVertex[i3], rVertex[i4]
                v1, v2, v3, v4 = vVertex[i1], vVertex[i2], vVertex[i3], vVertex[i4]

                n_vec = cross(r3 - r1, r2 - r4)
                normalUnitVec = n_vec / sqrt(dot(n_vec, n_vec))
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
    index = TEPanelIndex(model.boundMesh.sizes)
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
    mul!(Fa, model.segmentProps.aic3Ring, Γp)
    mul!(Fa, model.segmentProps.aic3Wake, Γw, 1.0, 1.0)
    
    nss = model.segmentProps.nSpanSegments
    sizes = model.boundMesh.sizes
    for (s, nc, ns) in sizes
        @batch for i in 1:nc
            for j in 1:ns
                i1, i2 = VertexIndex(s, i, j, sizes), VertexIndex(s, i+1, j, sizes)
                i3, i4 = VertexIndex(s, i, j+1, sizes), VertexIndex(s, i+1, j+1, sizes)
                v1, v2, v3, v4 = vVertex[i1], vVertex[i2], vVertex[i3], vVertex[i4]
                v = 0.375*(v1 + v2) + 0.125*(v3 + v4)
                m = SpanSegmentIndex(s, i, j, sizes)
                Fa[m] = ρ * Γs[m] * cross(Fa[m] + v, model.segmentProps.r[m])
            end
        end
        for j in 1:ns 
            m = SpanSegmentIndex(s, nc+1, j, sizes)
            Fa[m] = zero(eltype(Fa))
        end
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

function AeroSolve(vb::AbstractArray, model::AeroModel;
        ωb=SA[0.0, 0.0, 0.0],
        δc=zeros(eltype(vb), length(model.controlSurfaces)),
        dδc=zeros(eltype(vb), length(model.controlSurfaces)),
        rs=fill(zero(Point3{eltype(vb)}), model.boundMesh.sizes.totalVertices),
        vs=fill(zero(Point3{eltype(vb)}), model.boundMesh.sizes.totalVertices),
        ρ=1.225)

    T = promote_type(eltype(vb), eltype(ωb), eltype(δc), eltype(dδc))
    cache = CreateCacheArrays(model, T)
    AeroSolve!(cache, vb, ωb, δc, dδc, rs, vs, model, ρ)
    return cache
end

function AeroSolve(vb::AbstractArray{A}, ωb::AbstractArray{B},
        δc::AbstractArray{C}, model::AeroModel, ρ=1.225) where {A, B, C}
    T = promote_type(A, B, C)
    nVert = model.boundMesh.sizes.totalVertices
    rs = fill(zero(Point3{T}), nVert)
    vs = fill(zero(Point3{T}), nVert)
    dδc = zeros(T, length(model.controlSurfaces))
    cache = CreateCacheArrays(model, T)
    AeroSolve!(cache, vb, ωb, δc, dδc, rs, vs, model, ρ)
    return cache
end

function AeroSolve(vb::AbstractArray{A}, ωb::AbstractArray{B},
    δc::AbstractArray{C}, dδc::AbstractArray{D},
    rs::AbstractArray, vs::AbstractArray, model, ρ=1.225) where {A, B, C, D}

    AeroSolve(vb, model; ωb=ωb, δc=δc, dδc=dδc, rs=rs, vs=vs, ρ=ρ)
end

function AeroSolve!(cache, vb, ωb, δc, dδc, rs, vs, model, ρ=1.225)
    AddSteadyKinematics!(cache.rVertex, cache.vVertex, vb, ωb, δc, dδc, rs, vs, model)
    CalculateNormalwash!(cache.b, cache.rVertex, cache.vVertex, model)
    SolveCirculation!(cache.Γp, cache.Γw, cache.Γs, cache.b, model)
    CalculateAerodynamicForce!(cache.Fa, cache.Γp, cache.Γw, cache.Γs, cache.vVertex, model, ρ)
    return nothing
end

function GetTotalForces(cache, model)
    rRef = model.modelProperties.rRef
    return IntegrateLoad(cache.Fa, cache.Ra, rRef)
end

function GetStabilityCoefficients(cache, vb, ρ, model)
    mProps = model.modelProperties
    F, M = GetTotalForces(cache, model)
    α, _, _ = AerodynamicAngles(vb)    
    Fstab = GeometryToStabilityAxis(F, α, mProps)
    Mstab = GeometryToStabilityAxis(M, α, mProps)
    CFstab, CMstab = ConvertToCoefficients(Fstab, Mstab, vb, ρ, mProps)
    return CFstab, CMstab
end

function MonitorPointLoads!(Fmp, cache, model::SteadyAeroModel2D)
    for (i, mp) in enumerate(model.monitorPoints)
        indices = mp.segmentIndices
        Fvec = @view cache.Fa[indices]
        Rvec = @view cache.Ra[indices]
        F, M = IntegrateLoad(Fvec, Rvec, mp.origin)
        Fmp[6(i-1)+1:6(i-1)+6] .= [mp.orientation * F; mp.orientation * M]
    end
    return nothing
end

function GetStabilityDerivatives(model::SteadyAeroModel2D)
    mProps = model.modelProperties
    b, c = mProps.b, mProps.c
    V, ρ = 1.0, 1.0 
    
    function eval_coeffs(x::AbstractVector{T}) where T
        α, β, p_n, q_n, r_n = x[1], x[2], x[3], x[4], x[5]
        δc = x[6:end]
        ωb = SVector{3, T}(p_n * 2V / b, q_n * 2V / c, r_n * 2V / b)
        vb = BodyVelocity(V, α, β)
        cache = AeroSolve(vb, ωb, δc, model, ρ)
        CFstab, CMstab = GetStabilityCoefficients(cache, vb, ρ, model)
        return SVector{6, T}(CFstab[1], CFstab[2], CFstab[3], CMstab[1], CMstab[2], CMstab[3])
    end
    
    x0 = zeros(Float64, 5 + length(model.controlSurfaces))
    cfg = ForwardDiff.JacobianConfig(eval_coeffs, x0, ForwardDiff.Chunk{2}())
    dC_dx = ForwardDiff.jacobian(eval_coeffs, x0, cfg)
    return dC_dx
end

# FMI Bridge Implementation for SteadyAeroModel2D
function AeroPanels.AllocateFMUCaches(model::SteadyAeroModel2D{M}, ::Type{T}) where {M, T}
    return CreateCacheArrays(model, T)
end

function AeroPanels.InitializeFMU!(array::AbstractFMUArray, cache, model::SteadyAeroModel2D, t::Float64; start_from_trim::Bool=false)
    return nothing
end

function AeroPanels.EvaluateDerivatives!(du::AbstractVector, array::AbstractFMUArray, cache, model::SteadyAeroModel2D, t::Float64)
    return nothing
end

function AeroPanels.EvaluateDerivatives!(array::AbstractFMUArray, cache, model::SteadyAeroModel2D, t::Float64)
    return nothing
end

function AeroPanels.EvaluateOutputs!(array::AbstractFMUArray, cache, model::SteadyAeroModel2D, t::Float64)
    vb = SVector{3}(array.vel)
    ωb = SVector{3}(array.omega)
    δc = array.deltaC
    dδc = zeros(eltype(δc), length(δc))
    
    nVert = model.boundMesh.sizes.totalVertices
    T_el = eltype(vb)
    rs = reinterpret(Point3{T_el}, array.structDisplacement)
    vs = reinterpret(Point3{T_el}, array.structVelocity)
    ρ = array.rho[1]
    
    AeroSolve!(cache, vb, ωb, δc, dδc, rs, vs, model, ρ)
    
    Fg, Mg = GetTotalForces(cache, model)
    mProps = model.modelProperties
    Fb = GeometryToBodyAxis(Fg, mProps)
    Mb = GeometryToBodyAxis(Mg, mProps)
    array.forces .= Fb
    array.moments .= Mb
    
    α, _, _ = AerodynamicAngles(vb)
    CFbody, CMbody = ConvertToCoefficients(Fb, Mb, vb, ρ, mProps)
    array.coeffsBody .= (CFbody..., CMbody...)
    
    Fstab = GeometryToStabilityAxis(Fg, α, mProps)
    Mstab = GeometryToStabilityAxis(Mg, α, mProps)
    CFstab, CMstab = ConvertToCoefficients(Fstab, Mstab, vb, ρ, mProps)
    array.coeffsStab .= (CFstab..., CMstab...)
    
    if length(model.monitorPoints) > 0
        MonitorPointLoads!(array.monitorPointLoads, cache, model)
    end
    
    for (i, fvec) in enumerate(cache.Fa)
        array.nodalForces[3*(i-1)+1] = fvec[1]
        array.nodalForces[3*(i-1)+2] = fvec[2]
        array.nodalForces[3*(i-1)+3] = fvec[3]
    end
    
    return nothing
end
