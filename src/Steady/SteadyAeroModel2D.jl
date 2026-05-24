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
    controlSurfaces::Vector{ControlSurface{T}}
    monitorPoints::Vector{MonitorPoint{T}}
end

function AeroModel2D(surfaces::Vector{AeroSurface2D{T}}, props::AeroModelProperties{T};
    controlSurfaces=ControlSurface{T}[],
    monitorPoints=MonitorPoint{T}[]) where T
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
    return AeroModel2D(mesh, ringMesh, wakeMesh, sizes, wakeSizes, panelProperties,
     props, AIC, segmentProps, controlSurfaces, monitorPoints)
end

AerodynamicLoadLocation(model::AeroModel2D) = model.segmentProps.mid