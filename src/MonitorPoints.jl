struct MonitorPoint{T}
    name::String
    origin::SVector{3, T}
    orientation::SMatrix{3, 3, T, 9} 
    panelIndices::Vector{Int}       
    segmentIndices::Vector{Int}     
end

function CreateMonitorPoint(name, origin::SVector{3, T}, s::Int, nc::Tuple{Int,Int}, ns::Tuple{Int,Int},surfaces;
        orientation::SMatrix{3, 3, T, 9}=SMatrix{3,3}(1.,0,0, 0,1,0, 0, 0,1)) where T
    
    sizes = Sizes([size(surface) for surface in surfaces])
    nss = sizes.totalSpanSegments
    nts = sizes.totalSpanSegments + sizes.totalChordSegments
    spanIndices = [SpanSegmentIndex(s, i, j, sizes) for i in nc[1]:nc[2] for j in ns[1]:ns[2]]
    chordIndices = [nss + ChordSegmentIndex(s, i, j, sizes) for i in nc[1]:nc[2] for j in ns[1]:ns[2]]
    panelIndices = [nts + PanelIndex(s, i, j, sizes) for i in nc[1]:(nc[2]-1) for j in ns[1]:(ns[2]-1)]
    return MonitorPoint(name, origin, orientation, panelIndices, [spanIndices;chordIndices])
end

function IntegrateLoad(Fvec::AbstractArray{<:AbstractArray{A}},
                       Rvec::AbstractArray{<:AbstractArray{B}},
                       Rref::AbstractArray{C}) where {A, B, C}
    n = length(Fvec)
    n == length(Rvec) || throw(DimensionMismatch("Fvec and Rvec must have the same length ($n != $(length(Rvec)))"))

    F_total = zero(SVector{3, A})
    M_total = zero(SVector{3, promote_type(A, B, C)})

    @batch reduction=((+, F_total), (+, M_total)) for i in eachindex(Fvec)
        F_total += Fvec[i]
        M_total += cross(Rvec[i] - Rref, Fvec[i])
    end
    return F_total, M_total
end

function MonitorPointLoads(cache, model)
    Fmp = zeros(eltype(cache.b), 6 * length(model.monitorPoints))
    MonitorPointLoads!(Fmp, cache, model)
    return Fmp
end

struct MPDefinition{T}
    name::String
    surface::Int
    origin::SVector{3, T}
    nc::Tuple{Int, Int}
    ns::Tuple{Int, Int} 
    orientation::SMatrix{3, 3, T, 9}  
end
function MPDefinition(name, surface, origin, nc, ns, orientation=SMatrix{3,3}(1.,0,0, 0,1,0, 0, 0,1))
    MPDefinition(name, surface, origin, nc, ns, orientation)
end

function MonitorPoint(mp::MPDefinition, sizes)
    nss = sizes.totalSpanSegments
    nts = sizes.totalSpanSegments + sizes.totalChordSegments
    nc,ns = mp.nc, mp.ns

    spanIndices = [SpanSegmentIndex(mp.surface, i, j, sizes) for i in nc[1]:nc[2] for j in ns[1]:(ns[2]-1)]
    chordIndices = [nss + ChordSegmentIndex(mp.surface, i, j, sizes) for i in nc[1]:(nc[2]-1) for j in ns[1]:ns[2]]
    panelIndices = [nts + PanelIndex(mp.surface, i, j, sizes) for i in nc[1]:(nc[2]-1) for j in ns[1]:(ns[2]-1)]
    return MonitorPoint(mp.name, mp.origin, mp.orientation, panelIndices, [spanIndices;chordIndices])
end