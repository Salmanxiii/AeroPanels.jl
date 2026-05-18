struct MonitorPoint{T}
    name::String
    origin::SVector{3, T}
    orientation::SMatrix{3, 3, T, 9} 
    panelIndices::Vector{Int}       
    segmentIndices::Vector{Int}     
end

function MonitorPointLoads(F_vec, R_vec, mp::MonitorPoint)
    F_global = sum(F_vec[i] for i in indices)
    M_global = sum(cross(R_vec[i] - mp.origin, F_vec[i]) for i in indices) 
    F_local = mp.orientation * F_global
    M_local = mp.orientation * M_global
    return [F_local; M_local]
end