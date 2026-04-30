abstract type AeroModel end


################################## Aero Modeling #########################################

"""
$(TYPEDEF)

A 2D aerodynamic model representing a collection of lifting surfaces and their steady wake.

$(TYPEDFIELDS)
"""
struct AeroModel2D{T<:Real} <: AeroModel
    mesh::GeometryBasics.Mesh{3, T}
    ringMesh::GeometryBasics.Mesh{3, T}
    wakeMesh::GeometryBasics.Mesh{3, T}
    sizes::Sizes
    wakeSizes::Sizes
    panelProperties::PanelProperties{T}
    modelProperties::AeroModelProperties{T}
    AIC::LU{T, Matrix{T}, Vector{Int64}}
    segmentProps::SegmentProperties{T}
end

function AeroModel2D(surfaces::Vector{AeroSurface2D{T}}, props::AeroModelProperties{T}) where T
    mesh, sizes = CreateAeroMesh(surfaces)
    ringMesh = RingMesh(mesh, sizes)
    wakeMesh, wakeSizes = FlatWakeMesh(ringMesh, sizes, props)
    panelProperties = PanelProperties(sizes.totalPanels, mesh, FlowAxis(props))

    # AIC Calculation
    AIC = Influence(panelProperties.rCollocation, panelProperties.normal, ringMesh, props.symmXZ);
    AICwake = Influence(panelProperties.rCollocation, panelProperties.normal, wakeMesh, props.symmXZ);
    TEindices = TEPanelIndex(sizes)
    # Add Wake Influence
    AIC[:, TEindices] .+= AICwake
    AIC = lu(AIC)

    segmentProps = ProcessSegments(ringMesh, sizes, wakeMesh, wakeSizes, props.symmXZ)
    return AeroModel2D(mesh, ringMesh, wakeMesh, sizes, wakeSizes, panelProperties, props, AIC, segmentProps)
end

################################## Solution #########################################

"""
$(TYPEDEF)

The solution to a steady aerodynamic problem.

$(TYPEDFIELDS)
"""
struct SteadySolution{T}
    forceBody::Point3{T}
    forceStability::Point3{T}
    forceUnsteady::Point3{T}
    coefficientBody::Point3{T}
    coefficientStability::Point3{T}
    forceVecSeg::Vector{Point3{T}} # Force on segments
    CL::T
    CD::T
end

function SteadySolution(Fs, vb, aeroModel, ρ=1.225, Fu=[Point3(0.,0.,0.)])
    S = aeroModel.modelProperties.S
    Fu_sum = sum(Fu)
    F = sum(Fs) + Fu_sum
    Fbody = GeometryToBodyAxis(F, aeroModel.modelProperties)
    α, β, V = AerodynamicAngles(vb)
    Fstability = BodyFixedToStabilityAxis(Fbody, α)
    Q = 0.5 * dot(vb,vb) * ρ
    CL = -Fstability[3] / Q / S
    CD = -Fstability[1] / Q / S

    T = typeof(CL)
    return SteadySolution{T}(Fbody, Fstability, Fu_sum, Fbody/Q/S, Fstability/Q/S, Fs, CL, CD)
end

# 1. Compact Mode (for arrays/interpolation)
function Base.show(io::IO, sol::SteadySolution{T}) where T
    c = sol.coefficientStability
    print(io, "SteadySolution{$T}(CD=$(round(sol.CD, digits=4)), CY=$(round(c[2], digits=4)), CL=$(round(sol.CL, digits=4)))")
end

################################## Solution #########################################
"""
$(SIGNATURES)

Calculate the normal wash on each panel in-place into vector `b`.
"""
function NormalWash!(b::AbstractVector, vb::AbstractVector, model::AeroModel)
    @batch for i in 1:length(b)
        b[i] = -dot(vb, model.panelProperties.normal[i])
    end
    return nothing
end

"""
$(SIGNATURES)

Calculate and return the normal wash vector `b`.
"""
function NormalWash(vb::AbstractVector, model::AeroModel)
    b = Vector{eltype(vb)}(undef, model.sizes.totalPanels)
    NormalWash!(b, vb, model)
    return b
end

"""
$(SIGNATURES)

Solve for panel, wake, and segment circulations in-place.
"""
function Circulation!(Γp::AbstractVector, Γw::AbstractVector, Γs::AbstractVector, b::AbstractVector, model::AeroModel)
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

Solve for and return (Γp, Γw, Γs).
"""
function Circulation(b::AbstractVector, model::AeroModel)
    T = eltype(b)
    Γp = Vector{T}(undef, model.sizes.totalPanels)
    Γw = Vector{T}(undef, model.wakeSizes.totalPanels)
    Γs = Vector{T}(undef, model.segmentProps.nTotalSegments)
    Circulation!(Γp, Γw, Γs, b, model)
    return Γp, Γw, Γs
end

"""
$(SIGNATURES)

Calculate aerodynamic forces on segments in-place into vector `Fa`.
"""
function AerodynamicForce!(Fa::AbstractVector, Γp::AbstractVector, Γw::AbstractVector, Γs::AbstractVector, vb::AbstractVector, model::AeroModel; ρ=1.225)
    # Calculate induced velocity at segments
    # Non-allocating version of equation v = AIC3r * Γr + AIC3w * Γw
    Vi = model.segmentProps.aic3Ring * Γp
    mul!(Vi, model.segmentProps.aic3Wake, Γw, 1.0, 1.0)
    
    @batch for i in 1:length(Fa)
        Fa[i] = ρ * Γs[i] * cross(Vi[i] + vb, model.segmentProps.r[i])
    end
    return nothing
end

"""
$(SIGNATURES)

Calculate and return the aerodynamic force vector `Fa` on segments.
"""
function AerodynamicForce(Γp::AbstractVector, Γw::AbstractVector, Γs::AbstractVector, vb::AbstractVector, model::AeroModel; ρ=1.225)
    T = eltype(Γp)
    Fa = Vector{Point3{T}}(undef, model.segmentProps.nTotalSegments)
    AerodynamicForce!(Fa, Γp, Γw, Γs, vb, model, ρ=ρ)
    return Fa
end

"""
$(SIGNATURES)

Solve the steady aerodynamic problem for a given body velocity `vb`.
Returns a [`SteadySolution`](@ref).
"""
function AeroSolve(vb::AbstractVector, model::AeroModel2D, ρ = 1.225)
    b = NormalWash(vb, model)
    Γp, Γw, Γs = Circulation(b, model)
    Fa = AerodynamicForce(Γp, Γw, Γs, vb, model, ρ=1.225)
    return SteadySolution(Fa, vb, model, ρ)
end

"""
$(SIGNATURES)

Solve the steady aerodynamic problem for a given speed `V` and angle of attack `α` (degrees).
"""
function AeroSolve(V, α::Real, model::AeroModel2D, ρ = 1.225)
    vb = BodyVelocity(V, deg2rad(α))
    return AeroSolve(vb, model, ρ)
end

AerodynamicLoadLocation(model::AeroModel2D) = model.segmentProps.mid