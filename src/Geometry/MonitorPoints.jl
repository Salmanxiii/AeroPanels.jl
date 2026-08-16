struct MonitorPoint{T}
    name::String
    origin::SVector{3, T}
    orientation::SMatrix{3, 3, T, 9} 
    panelIndices::Vector{Int}       
    segmentIndices::Vector{Int}     
end

function CreateMonitorPoint(name::String, origin::SVector{3, T}, panelIndices::Vector{Int}, segmentIndices::Vector{Int};
        orientation::SMatrix{3, 3, T, 9}=SMatrix{3,3,T}(1.,0,0, 0,1,0, 0, 0,1)) where T
    return MonitorPoint{T}(name, origin, orientation, panelIndices, segmentIndices)
end

function CreateMonitorPoint(mesh::ThinAeroMesh{T}, name::String; 
                            surfaceIdx::Int, nc::Tuple{Int,Int}, ns::Tuple{Int,Int},
                            origin::Union{Nothing, SVector{3, T}}=nothing,
                            orientation::SMatrix{3, 3, T, 9}=SMatrix{3,3,T}(1.,0,0, 0,1,0, 0, 0,1)) where T
    sizes = mesh.sizes
    nss = sizes.totalSpanSegments
    nts = sizes.totalSpanSegments + sizes.totalChordSegments
    
    spanIndices = [SpanSegmentIndex(surfaceIdx, i, j, sizes) for i in nc[1]:nc[2] for j in ns[1]:(ns[2]-1)]
    chordIndices = [nss + ChordSegmentIndex(surfaceIdx, i, j, sizes) for i in nc[1]:(nc[2]-1) for j in ns[1]:ns[2]]
    panelIndices = [nts + PanelIndex(surfaceIdx, i, j, sizes) for i in nc[1]:(nc[2]-1) for j in ns[1]:(ns[2]-1)]
    
    segmentIndices = [spanIndices; chordIndices]
    
    calcOrigin = isnothing(origin) ? mesh.rCollocation[panelIndices[1]] : origin
    
    return MonitorPoint{T}(name, calcOrigin, orientation, panelIndices, segmentIndices)
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
