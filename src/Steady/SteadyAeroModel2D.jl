################################## Aero Modeling #########################################

"""
$(TYPEDEF)

A 2D aerodynamic model representing a collection of lifting surfaces and their steady wake.

$(TYPEDFIELDS)
"""
struct AeroModel2D{T<:Real} <: AeroModel
    mesh::Mesh{3, T, QuadFace{Int64}, (:position,), Tuple{Vector{Point{3, T}}}, Vector{QuadFace{Int64}}}
    ringMesh::Mesh{3, T, QuadFace{Int64}, (:position,), Tuple{Vector{Point{3, T}}}, Vector{QuadFace{Int64}}}
    wakeMesh::Mesh{3, T, QuadFace{Int64}, (:position,), Tuple{Vector{Point{3, T}}}, Vector{QuadFace{Int64}}}
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
    controlSurfaces=CSDefinition[],
    monitorPoints=MPDefinition[]) where T
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
    controlSurfaces2 = [ControlSurface(cs, sizes, mesh) for cs in controlSurfaces]
    monitorPoints2 = [MonitorPoint(mp, sizes) for mp in monitorPoints]
    
    return AeroModel2D(mesh, ringMesh, wakeMesh, sizes, wakeSizes, panelProperties,
     props, AIC, segmentProps, controlSurfaces2, monitorPoints2)
end

AerodynamicLoadLocation(model::AeroModel2D) = model.segmentProps.mid