abstract type AeroModel end


################################## Aero Modeling #########################################

struct AeroModel2D{T<:Real} <: AeroModel
    mesh::GeometryBasics.Mesh{3, T}
    ringMesh::GeometryBasics.Mesh{3, T}
    wakeMesh::GeometryBasics.Mesh{3, T}
    sizes::Sizes
    wakeSizes::Sizes
    panelProperties::PanelProperties{T}
    modelProperties::AeroModelProperties{T}
    AIC::LU{T, Matrix{T}, Vector{Int64}}
    segmentSpan::SegmentProperties{T}
    segmentChord::SegmentProperties{T}
end

function AeroModel2D(surfaces::Vector{AeroSurface2D{T}}, props::AeroModelProperties{T}) where T
    mesh, sizes = CreateAeroMesh(surfaces)
    ringMesh = RingMesh(mesh, sizes)
    wakeMesh, wakeSizes = FlatWakeMesh(ringMesh, sizes, props)
    panelProperties = PanelProperties(sizes.totalPanels, mesh, FlowAxis(props))

    AIC = SteadyWakeInfluence(panelProperties.rCollocation, panelProperties.normal,
     ringMesh, wakeMesh, sizes, props.symmXZ)
    AIC = lu(AIC)
    segmentSpan, segmentChord = ProcessSegments(ringMesh, sizes, wakeMesh, wakeSizes, props.symmXZ)

    return AeroModel2D(mesh, ringMesh, wakeMesh, sizes, wakeSizes, panelProperties, props, AIC, segmentSpan, segmentChord)
end

################################## Solution #########################################

struct SteadySolution{T}
    forceBody::Point3{T}
    forceStability::Point3{T}
    forceUnsteady::Point3{T}
    coefficientBody::Point3{T}
    coefficientStability::Point3{T}
    forceVecSpan::Vector{Point3{T}} # Force on spanwise segments
    forceVecChord::Vector{Point3{T}}
    CL::T
    CD::T
end

function SteadySolution(Fs, Fc, vb, aeroModel, ρ=1.225, Fu=[Point3(0.,0.,0.)])
    S = aeroModel.modelProperties.S
    Fu_sum = sum(Fu)
    F = sum(Fs) + sum(Fc) + Fu_sum
    Fbody = GeometryToBodyAxis(F, aeroModel.modelProperties)
    α, β, V = AerodynamicAngles(vb)
    Fstability = BodyFixedToStabilityAxis(Fbody, α)
    Q = 0.5 * dot(vb,vb) * ρ
    CL = -Fstability[3] / Q / S
    CD = -Fstability[1] / Q / S

    T = typeof(CL)
    return SteadySolution{T}(Fbody, Fstability, Fu_sum, Fbody/Q/S, Fstability/Q/S, Fs, Fc, CL, CD)
end

# 1. Compact Mode (for arrays/interpolation)
function Base.show(io::IO, sol::SteadySolution{T}) where T
    c = sol.coefficientStability
    print(io, "SteadySolution{$T}(CD=$(round(sol.CD, digits=4)), CY=$(round(c[2], digits=4)), CL=$(round(sol.CL, digits=4)))")
end

################################## Solution #########################################

function AICSolve(Vvec::Vector, aeroModel::AeroModel2D)
    RHS = -[dot(v, nrml) for (v, nrml) in zip(Vvec, aeroModel.panelProperties.normal)]
    Γp = aeroModel.AIC\RHS
    index = TEPanelIndex(aeroModel.sizes)
    Γw = @view Γp[index]
    return Γp, Γw
end

function AerodynamicForces(Γp, Γw, aeroModel::AeroModel2D, vb, ρ)
    # Γs, Γc = SegmentCirculation(Γp, aeroModel.sizes);
    Fs = SegmentForce(Γp, Γw, vb, ρ, aeroModel.segmentSpan)
    Fc = SegmentForce(Γp, Γw, vb, ρ, aeroModel.segmentChord)
    return Fs, Fc
end

function AeroSolve(V, vb::AbstractVector, aeroModel::AeroModel2D, ρ = 1.225)

    velVec = [vb for i in 1:aeroModel.sizes.totalPanels]
    Γp, Γw = AICSolve(velVec, aeroModel)
    Fs, Fc = AerodynamicForces(Γp, Γw, aeroModel, vb, ρ)

    return SteadySolution(Fs, Fc, vb, aeroModel, ρ)
end

function AeroSolve(V, α::Real, aeroModel::AeroModel2D, ρ = 1.225)
    vb = BodyVelocity(V, deg2rad(α))
    return AeroSolve(V, vb, aeroModel, ρ)
end

