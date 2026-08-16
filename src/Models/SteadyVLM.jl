################################## Aero Modeling #########################################

"""
$(TYPEDEF)

A 2D aerodynamic model representing a collection of lifting surfaces and their steady wake.

$(TYPEDFIELDS)
"""
@fmu_model struct SteadyAeroModel2D{T<:Real} <: AbstractFMUModel{T}
    boundMesh::ThinAeroMesh{T}
    wakeMesh::FixedWakeMesh{T}
    modelProperties::AeroModelProperties{T}
    AIC::LU{T, Matrix{T}, Vector{Int64}}
    segmentProps::SegmentProperties{T}
    controlSurfaces::Vector{ControlSurface{T}}
    monitorPoints::Vector{MonitorPoint{T}}
    
    [States]
    # No states for steady model

    [Derivatives]
    # No derivatives for steady model

    [Inputs]
    pos                = 3
    vel                = 3
    omega              = 3
    euler              = 3
    ab                 = 3
    domega             = 3
    deltaC             = (m) -> length(m.controlSurfaces)
    structDisplacement = (m) -> 3 * m.boundMesh.sizes.totalVertices
    structVelocity     = (m) -> 3 * m.boundMesh.sizes.totalVertices
    structAcceleration = (m) -> 3 * m.boundMesh.sizes.totalVertices
    rho                = 1

    [Outputs]
    forces            = 3
    moments           = 3
    coeffsBody        = 6
    coeffsStab        = 6
    monitorPointLoads = (m) -> 6 * length(m.monitorPoints)
    nodalForces       = (m) -> 3 * TotalSegments(m.boundMesh)
end

function SteadyAeroModel2D(boundMesh::ThinAeroMesh{T}, wakeMesh::FixedWakeMesh{T}, props::AeroModelProperties{T};
    controlSurfaces=ControlSurface{T}[],
    monitorPoints=MonitorPoint{T}[]) where T
    
    rCollocation = GetCollocationPoints(boundMesh)
    normals = GetNormalVectors(boundMesh)
    
    AIC = Influence(rCollocation, normals, boundMesh.ringMesh, props.symmXZ);
    AICwake = Influence(rCollocation, normals, wakeMesh.mesh, props.symmXZ);
    TEindices = TEPanelIndex(boundMesh.sizes)
    
    # Add Wake Influence
    AIC[:, TEindices] .+= AICwake
    AIC = lu(AIC)

    segmentProps = ProcessSegments(boundMesh.ringMesh, boundMesh.sizes, wakeMesh.mesh, wakeMesh.wakeSizes, props.symmXZ)
    
    return SteadyAeroModel2D(boundMesh, wakeMesh, props, AIC, segmentProps, controlSurfaces, monitorPoints)
end
