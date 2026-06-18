struct ControlSurface{T}
    name::String
    origin::SVector{3, T}
    hinge::SVector{3, T}
    index::Vector{Int}         # Vertex indices (for physical mesh displacement)
    panelIndices::Vector{Int}  # Panel indices (for aerodynamic boundary conditions)
end

function CreateControlSurface(name, origin, hingeAxis, s::Int, nc::Tuple{Int,Int}, ns::Tuple{Int,Int}, surfaces)
    sizes = Sizes([size(surface) for surface in surfaces])
    indices = [VertexIndex(s, i, j, sizes) for i in nc[1]:nc[2] for j in ns[1]:ns[2]]
    pIndices = [PanelIndex(s, i, j, sizes) for i in nc[1]:(nc[2]-1) for j in ns[1]:(ns[2]-1)]
    return ControlSurface(name, origin, normalize(hingeAxis), indices, pIndices)
end

function CSDisplacement!(vertices::Vector, δ, cs::ControlSurface)
    cδ = cos(δ)
    sδ = sin(δ)
    @batch for i in cs.index
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
    @batch for i in cs.index
        r = vertices[i] - cs.origin
        velocities[i] = -dδ * cross(cs.hinge, r)
    end
end

function CSVelocity!(velocities::Vector, dδc, model::AeroModel)
    vertices = coordinates(model.mesh)
    for (cs, dδ) in zip(model.controlSurfaces, dδc)
        CSVelocity!(velocities, dδ, vertices, cs)
    end
end

function CSAcceleration!(accelerations::Vector, ddδc, model::AeroModel)
    vertices = coordinates(model.mesh)
    for (cs, ddδ) in zip(model.controlSurfaces, ddδc)
        CSVelocity!(accelerations, ddδ, vertices, cs)
    end
end

struct CSDefinition
    name::String
    surface::Int
    nc::Tuple{Int, Int}
    ns::Tuple{Int, Int}
end

function ControlSurface(cs::CSDefinition, sizes::Sizes, mesh)
    r = coordinates(mesh)
    i1 = VertexIndex(cs.surface, cs.nc[1], cs.ns[1], sizes)
    i2 = VertexIndex(cs.surface, cs.nc[1], cs.ns[2], sizes)
    hingeAxis = r[i2] - r[i1]

    indices = [VertexIndex(cs.surface, i, j, sizes) for i in cs.nc[1]:cs.nc[2] for j in cs.ns[1]:cs.ns[2]]
    pIndices = [PanelIndex(cs.surface, i, j, sizes) for i in cs.nc[1]:(cs.nc[2]-1) for j in cs.ns[1]:(cs.ns[2]-1)]
    return ControlSurface(cs.name, SVector(r[i1]), SVector(normalize(hingeAxis)), indices, pIndices)
end