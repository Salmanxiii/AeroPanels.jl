############################# Model Properties ######################################
"""
$(TYPEDEF)

Aerodynamic model properties.

$(TYPEDFIELDS)
"""
@kwdef struct AeroModelProperties{T}
    c::T
    b::T
    S::T
    symmXZ::Bool = false
    symmXY::Bool = false
#    flowAxis::SVector{3,T} = SA[1.0, 0., 0.]
    bodyFixedCS::SMatrix{3,3,T,9} = SMatrix{3,3}(-1.,0,0, 0,1,0, 0, 0,-1)
    rRef::SVector{3, T} = SVector{3}(0., 0., 0.)
end

FlowAxis(s::AeroModelProperties) = -s.bodyFixedCS[:, 1]

############################# ThinAeroMesh ######################################
"""
$(TYPEDEF)

Properties of the panels in the aerodynamic mesh.

$(TYPEDFIELDS)
"""
struct ThinAeroMesh{T} <: AbstractAeroMesh{T}
    surfaces::Vector{AeroSurface2D{T}}
    mesh::GeometryBasics.Mesh{3, T, QuadFace{Int64}}
    ringMesh::GeometryBasics.Mesh{3, T, QuadFace{Int64}}
    sizes::Sizes
    normals::Vector{Point3{T}}
    rMid::Vector{Point3{T}}
    rCollocation::Vector{Point3{T}}
    rControl::Vector{Point3{T}}
    areas::Vector{T}
    spans::Vector{T}
    symmetryXZ::Bool
end

function ThinAeroMesh(surfaces::Vector{AeroSurface2D{T}}; symmetryXZ::Bool=false, flowAxis=SA[1.0, 0.0, 0.0]) where T
    mesh, sizes = CreateAeroMesh(surfaces)
    ringMesh = RingMesh(mesh, sizes)
    
    totalPanels = sizes.totalPanels
    normals = Vector{Point3{T}}(undef, totalPanels)
    rMid = Vector{Point3{T}}(undef, totalPanels)
    rCollocation = Vector{Point3{T}}(undef, totalPanels)
    rControl = Vector{Point3{T}}(undef, totalPanels)
    areas = Vector{T}(undef, totalPanels)
    spans = Vector{T}(undef, totalPanels)
    
    # Multithreaded for loop to calculate the properties of each panel
    @batch for i in 1:length(mesh)
        r1, r2, r3, r4 = mesh[i]
        nVec = cross(r3 - r1, r2 - r4)
        normals[i] = normalize(nVec)
        area = norm(nVec) / 2
        areas[i] = area
        LEmid = (r1 + r2) / 2
        TEmid = (r3 + r4) / 2
        chordMid = dot(flowAxis, TEmid - LEmid)
        spans[i] = area / chordMid
        rControl[i] = 0.75 * LEmid + 0.25 * TEmid
        rMid[i] = (LEmid + TEmid) / 2
        rCollocation[i] = 0.25 * LEmid + 0.75 * TEmid
    end
    
    return ThinAeroMesh{T}(surfaces, mesh, ringMesh, sizes, normals, rMid, rCollocation, rControl, areas, spans, symmetryXZ)
end

GetCollocationPoints(aeroMesh::ThinAeroMesh) = aeroMesh.rCollocation
GetMidpoints(aeroMesh::ThinAeroMesh)    = aeroMesh.rMid
GetNormalVectors(aeroMesh::ThinAeroMesh)     = aeroMesh.normals
GetPanelAreas(aeroMesh::ThinAeroMesh)        = aeroMesh.areas
GetBoundPanelCount(aeroMesh::ThinAeroMesh)   = aeroMesh.sizes.totalPanels
GetBoundMesh(aeroMesh::ThinAeroMesh)         = aeroMesh.mesh

TotalSegments(aeroMesh::ThinAeroMesh) = aeroMesh.sizes.totalSpanSegments + aeroMesh.sizes.totalChordSegments
TotalPanels(aeroMesh::ThinAeroMesh)   = aeroMesh.sizes.totalPanels

struct FixedWakeMesh{T} <: AbstractAeroMesh{T}
    mesh::GeometryBasics.Mesh{3, T, QuadFace{Int64}}
    wakeSizes::Sizes
    wakeSpacing::Vector{T}
end

GetWakeMesh(aeroMesh::FixedWakeMesh) = aeroMesh.mesh

function FixedWakeMesh(ringMesh::GeometryBasics.Mesh{3, T}, sizes, props::AeroModelProperties{T}; nWake::Int64=1, wakeLength::T=20., wakeExpansionR::Float64=1.0) where T
    wakeSizes = Sizes([(nWake, ns) for (s, nc, ns) in sizes])
    
    L_wake = wakeLength * props.c
    wakeSpacing = zeros(T, nWake)
    
    if wakeExpansionR ≈ 1.0
        fill!(wakeSpacing, L_wake / nWake)
    else
        dx_first = 2.0 * L_wake / (nWake * (1.0 + wakeExpansionR))
        d = dx_first * (wakeExpansionR - 1.0) / (nWake - 1.0)
        for k in 1:nWake
            wakeSpacing[k] = dx_first + (k - 1) * d
        end
    end
    
    wakeOffsets = zeros(T, nWake + 1)
    for i in 2:(nWake + 1)
        wakeOffsets[i] = wakeOffsets[i - 1] + wakeSpacing[i - 1]
    end

    # Initialize Arrays
    wakeVertices = Vector{Point3{T}}(undef, wakeSizes.totalVertices)
    wakeFaces = Vector{GeometryBasics.QuadFace{Int}}(undef, wakeSizes.totalPanels)
    ringVertices = coordinates(ringMesh)
    flowDir = FlowAxis(props)

    for (s, nc, ns) in sizes
        @batch for i in 1:nWake+1
            rw = wakeOffsets[i] * flowDir
            for j in 1:ns+1
                rRingTE = ringVertices[VertexIndex(s, nc+1, j, sizes)]
                index = VertexIndex(s, i, j, wakeSizes)
                wakeVertices[index] = rRingTE + rw
            end
        end

        @batch for i in 1:nWake
            for j in 1:ns
                index = PanelIndex(s, i, j, wakeSizes)
                wakeFaces[index] = GeometryBasics.QuadFace(
                    VertexIndex(s, i,   j,   wakeSizes),
                    VertexIndex(s, i,   j+1, wakeSizes),
                    VertexIndex(s, i+1, j+1, wakeSizes),
                    VertexIndex(s, i+1, j,   wakeSizes))
            end
        end
    end
    wakeMesh = type_stable_mesh((position=wakeVertices,), wakeFaces, UnitRange{UInt32}[])
    return FixedWakeMesh{T}(wakeMesh, wakeSizes, wakeSpacing)
end
