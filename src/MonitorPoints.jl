struct MonitorPoint{T}
    name::String
    origin::SVector{3, T}
    orientation::SMatrix{3, 3, T, 9} 
    panelIndices::Vector{Int}       
    segmentIndices::Vector{Int}     
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

function MonitorPointLoads(FVec, RVec, mp::MonitorPoint)
    F, M = IntegrateLoad(view(Fvec, mp.indecs), view(Rvec, mp.indices), mp.origin)
    F_local = mp.orientation * F
    M_local = mp.orientation * M
    return [F_local; M_local]
end