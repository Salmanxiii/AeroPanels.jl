struct ControlSurface{T}
    name::String
    origin::SVector{3, T}
    hinge::SVector{3, T}
    vertexIndices::Vector{Int}
    panelIndices::Vector{Int}
end

# Method 1: General Constructor (For Unstructured / Imported Meshes using global indices)
function CreateControlSurface(name::String, origin::SVector{3,T}, hingeAxis::SVector{3,T}, 
                              vertexIndices::Vector{Int}, panelIndices::Vector{Int}) where T
    return ControlSurface{T}(name, origin, normalize(hingeAxis), vertexIndices, panelIndices)
end

# Method 2: Structured Constructor (Specific to ThinAeroMesh using surface/chord/span indices)
function CreateControlSurface(mesh::ThinAeroMesh{T}, name::String; 
                              surfaceIdx::Int, nc::Tuple{Int, Int}, ns::Tuple{Int, Int}, 
                              origin::Union{Nothing, SVector{3, T}}=nothing, 
                              hingeAxis::Union{Nothing, SVector{3, T}}=nothing) where T
    sizes = mesh.sizes
    
    # Generate structural indices perfectly aligned with panel boundaries
    vertexIndices = [VertexIndex(surfaceIdx, i, j, sizes) for i in nc[1]:(nc[2]+1) for j in ns[1]:(ns[2]+1)]
    panelIndices = [PanelIndex(surfaceIdx, i, j, sizes) for i in nc[1]:nc[2] for j in ns[1]:ns[2]]
    
    # Calculate origin and hinge axis if not provided
    r = coordinates(mesh.mesh)
    i1 = VertexIndex(surfaceIdx, nc[1], ns[1], sizes)
    i2 = VertexIndex(surfaceIdx, nc[1], ns[2]+1, sizes)
    
    calcOrigin = isnothing(origin) ? r[i1] : origin
    calcHinge = isnothing(hingeAxis) ? (r[i2] - r[i1]) : hingeAxis
    
    return ControlSurface{T}(name, SVector{3,T}(calcOrigin), SVector{3,T}(normalize(calcHinge)), vertexIndices, panelIndices)
end

function CSDisplacement!(vertices::Vector, δ, cs::ControlSurface)
    cδ = cos(δ)
    sδ = sin(δ)
    @batch for i in cs.vertexIndices
        r = vertices[i] - cs.origin
        v = cross(cs.hinge, r)
        vertices[i] = vertices[i] + sδ * v + (1 - cδ) * cross(cs.hinge, v)
    end
end

function CSDisplacement!(vertices, δc, model::AeroModel)
    for (cs, δ) in zip(model.controlSurfaces, δc)
        CSDisplacement!(vertices, δ, cs)
    end
end

function CSVelocity!(velocities::Vector, dδ, vertices::Vector, cs::ControlSurface)
    @batch for i in cs.vertexIndices
        r = vertices[i] - cs.origin
        velocities[i] = -dδ * cross(cs.hinge, r)
    end
end

function CSVelocity!(velocities::Vector, dδc, model::AeroModel)
    vertices = coordinates(GetBoundMesh(model))
    for (cs, dδ) in zip(model.controlSurfaces, dδc)
        CSVelocity!(velocities, dδ, vertices, cs)
    end
end

function CSAcceleration!(accelerations::Vector, ddδc, model::AeroModel)
    vertices = coordinates(GetBoundMesh(model))
    for (cs, ddδ) in zip(model.controlSurfaces, ddδc)
        CSVelocity!(accelerations, ddδ, vertices, cs)
    end
end