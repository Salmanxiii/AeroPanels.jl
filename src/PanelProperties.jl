
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
end

FlowAxis(s::AeroModelProperties) = -s.bodyFixedCS[:, 1]
############################# PanelProperties ######################################
"""
$(TYPEDEF)

Properties of the panels in the aerodynamic mesh.

$(TYPEDFIELDS)
"""
struct PanelProperties{T}
    normal::Vector{Point3{T}}
    rMid::Vector{Point3{T}}
    rCollocation::Vector{Point3{T}}
    rControl::Vector{Point3{T}}
    area::Vector{T}
    span::Vector{T}
end

function PanelProperties(totalPanels::Int, mesh::GeometryBasics.Mesh{3, T}, flowAxis::SVector) where T
    p = PanelProperties{T}(totalPanels)
    PanelProperties!(p, mesh, flowAxis)
    return p
end

# Constructor to initialize fields based on totalPanels
function PanelProperties{T}(totalPanels::Int) where T
    normal = Vector{Point3{T}}(undef, totalPanels)
    rMid = Vector{Point3{T}}(undef, totalPanels)
    rCollocation = Vector{Point3{T}}(undef, totalPanels) # 75% chord
    rControl = Vector{Point3{T}}(undef, totalPanels) # 25% chord
    area = Vector{T}(undef, totalPanels)
    span = Vector{T}(undef, totalPanels)
    return PanelProperties(normal, rMid, rCollocation, rControl, area, span)
end

function PanelProperties!(p::PanelProperties{T}, mesh::GeometryBasics.Mesh, flowAxis::SVector{3, T}) where T
    # Multithreaded for loop to calculate the properties of each panel
    @batch for i in 1:length(mesh)
        r1, r2, r3, r4 = mesh[i]
        nVec = cross(r3 - r1, r2 - r4)
        p.normal[i] = normalize(nVec)
        area = norm(nVec) / 2
        p.area[i] = area
        LEmid = (r1 + r2) / 2
        TEmid = (r3 + r4) / 2
        chordMid = dot(flowAxis, TEmid - LEmid)
        p.span[i] = area / chordMid
        p.rControl[i] = 0.75 * LEmid + 0.25 * TEmid
        p.rMid[i] = (LEmid + TEmid) / 2
        p.rCollocation[i] = 0.25 * LEmid + 0.75 * TEmid
    end
end
