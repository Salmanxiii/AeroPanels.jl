struct ControlSurface{T}
name:: String
origin::SVector{3, T}
hinge::SVector{3, T}
index::Vector{Int}
end

function CreateControlSurface(name, origin, hingeAxis, s::Int, nc::Tuple{Int,Int}, ns::Tuple{Int,Int}, surfaces)
    sizes = Sizes([size(surface) for surface in surfaces])
    indices = [VertexIndex(s, i, j, sizes) for i in nc[1]:nc[2] for j in ns[1]:ns[2]]
    return ControlSurface(name, origin, normalize(hingeAxis)  , indices)
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
        velocities[i] = dδ * cross(cs.hinge, r)
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

function NormalWash(vb, ab, ωb, dωb, δc, dδc, ddδc, vd, ad, model)
    n = NumberOfStates(model)
    b = zeros(n)
    db = copy(b)
    vertices = copy(coordinates(model.mesh))
    velocities = zeros(eltype(vertices), length(vertices))
    accelerations = zeros(eltype(vertices), length(vertices))

    CSDisplacement!(vertices, δc, model)
    CSVelocity!(velocities, dδc, model)
    CSAcceleration!(accelerations, ddδc, model)

    bodyAxis = model.modelProperties.bodyFixedCS # body axis coordinate
    CG = model.modelProperties.CG
    # velocities in geometry frame
    (vb_g, ab_g, ωb_g, dωb_g) = (bodyAxis' * component for component in (vb, ab, ωb, dωb))

    for i in eachindex(vertices)
        r = vertices[i] - CG
        velocities[i] = velocities[i] + vb_g + vd[i] + cross(ωb_g, r)
        accelerations[i] = accelerations[i] + ab_g + ad[i] + cross(dωb_g, r) # centripetal acceleration ignored
    end

    # Normalwash calculation: calculate new normals and interpolate velocity at collocation point
    for (s, nc, ns) in model.sizes
        @batch for i in 1:nc
            for j in 1:ns
                p = PanelIndex(s, i, j, model.sizes)
                i1, i2, i3, i4 = PanelVertexIndices(s, i, j, model.sizes)

                r1, r2, r3, r4 = vertices[i1], vertices[i2], vertices[i3], vertices[i4]
                v1, v2, v3, v4 = velocities[i1], velocities[i2], velocities[i3], velocities[i4]
                a1, a2, a3, a4 = accelerations[i1], accelerations[i2], accelerations[i3], accelerations[i4]

                normalVec = cross(r3 - r1, r2 - r4)
                normalNorm = norm(normalVec)
                normalUnitVec = normalVec / normalNorm
                v = 0.125*(v1 + v2) + 0.375*(v3 + v4)
                a = 0.125*(a1 + a2) + 0.375*(a3 + a4)
                b[p] = - dot(v, normalUnitVec)

                dNormalVec = cross(v3 - v1, r2 - r4) + cross(r3 - r1, v2 - v4)
                dNormalUnitVec = (dNormalVec - normalUnitVec * dot(normalUnitVec, dNormalVec)) / normalNorm
                db[p] = -(dot(a, normalUnitVec) + dot(v, dNormalUnitVec))
            end
        end
    end
    return b, db
end