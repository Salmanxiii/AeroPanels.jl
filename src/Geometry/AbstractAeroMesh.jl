abstract type AbstractAeroMesh{T} end

GetCollocationPoints(aeroMesh::AbstractAeroMesh) = error("GetCollocationPoints not implemented for $(typeof(aeroMesh))")
GetMidpoints(aeroMesh::AbstractAeroMesh)    = error("GetMidpoints not implemented for $(typeof(aeroMesh))")
GetNormalVectors(aeroMesh::AbstractAeroMesh)     = error("GetNormalVectors not implemented for $(typeof(aeroMesh))")
GetPanelAreas(aeroMesh::AbstractAeroMesh)        = error("GetPanelAreas not implemented for $(typeof(aeroMesh))")
GetBoundPanelCount(aeroMesh::AbstractAeroMesh)   = error("GetBoundPanelCount not implemented for $(typeof(aeroMesh))")
GetBoundMesh(aeroMesh::AbstractAeroMesh)         = error("GetBoundMesh not implemented for $(typeof(aeroMesh))")
GetWakeMesh(aeroMesh::AbstractAeroMesh)          = error("GetWakeMesh not implemented for $(typeof(aeroMesh))")
GetControlSurfaces(aeroMesh::AbstractAeroMesh)   = error("GetControlSurfaces not implemented for $(typeof(aeroMesh))")
GetMonitorPoints(aeroMesh::AbstractAeroMesh)     = error("GetMonitorPoints not implemented for $(typeof(aeroMesh))")
