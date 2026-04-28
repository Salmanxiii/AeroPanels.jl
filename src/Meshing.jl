abstract type WakeModel end

################################## Meshing #########################################

function CreateAeroMesh(surfaces::Vector{AeroSurface2D{T}}) where T
    sizes = Sizes([(surf.nc, surf.ns) for surf in surfaces])
    
    vertices = Vector{Point3{T}}(undef, sizes.totalVertices)
    faces = Vector{GeometryBasics.QuadFace{Int}}(undef, sizes.totalPanels)

    index1, index2 = 1, 1
    offset = 0
    for (s, nc, ns) in sizes
        X, Y, Z = surfaces[s].X, surfaces[s].Y, surfaces[s].Z
        
        for i in 1:nc+1, j in 1:ns+1
            vertices[index1] = Point3(X[i, j], Y[i, j], Z[i, j])
            index1 += 1
        end
        # Indices of each face/panel
        for i in 1:nc, j in 1:ns
            faces[index2] = GeometryBasics.QuadFace(
                offset + (i-1)*(ns+1)+j,
                offset + (i-1)*(ns+1)+j+1,
                offset + i*(ns+1)+j+1,
                offset + i*(ns+1)+j)
            index2 += 1
        end
        offset += (nc+1)*(ns+1)
    end
    mesh = GeometryBasics.Mesh(vec(vertices), vec(faces))
    return mesh, sizes
end

function RingMesh!(ringMesh, mesh, sizes)
    rVert = coordinates(mesh)
    rRing = coordinates(ringMesh)
    for (s, nc, ns) in sizes
        for i in 1:nc
            for j in 1:ns+1
                index = VertexIndex(s, i , j, sizes)
                index2 = VertexIndex(s, i+1 , j, sizes)
                rRing[index] = 0.75*rVert[index] + 0.25*rVert[index2]
            end
        end
        for j in 1:ns+1
            index = VertexIndex(s, nc , j, sizes)
            index2 = VertexIndex(s, nc+1 , j, sizes)
            r1, r2 = rVert[index], rVert[index2]
            rRing[index2] = r2 + 0.25*(r2 - r1)
        end
    end
end

function RingMesh(mesh, sizes)
    ringMesh = GeometryBasics.Mesh(similar(coordinates(mesh)), faces(mesh))
    RingMesh!(ringMesh, mesh, sizes)
    return ringMesh
end


################################## Wake Modeling #########################################

@kwdef struct SteadyWake <: WakeModel
    length::Real = 30.0
end

@kwdef struct FixedWake <: WakeModel
    n::Int
    length::Real = 30.0
end

function FlatWakeMesh(ringMesh::GeometryBasics.Mesh{3, T}, sizes, props::AeroModelProperties{T}; nWake::Int64=1, wakeLength::T=20.) where T
    wakeSizes = Sizes([(nWake, ns) for (s, nc, ns) in sizes])
    wakePanelLength = wakeLength*props.c / nWake
    # Initialize Arrays
    wakeVertices = Vector{Point3{T}}(undef, wakeSizes.totalVertices)
    wakeFaces = Vector{GeometryBasics.QuadFace{Int}}(undef, wakeSizes.totalPanels)
    ringVertices = coordinates(ringMesh)
    flowDir = FlowAxis(props)

    for (s, nc, ns) in sizes
        @batch for i in 1:nWake+1
            rw = (i-1)*wakePanelLength*flowDir
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
    return GeometryBasics.Mesh(wakeVertices, wakeFaces), wakeSizes
end
