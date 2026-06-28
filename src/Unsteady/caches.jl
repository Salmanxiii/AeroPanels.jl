
# --- Cache Constructors ---

"""
    InputCache(model::UnsteadyAeroModel2D{M}, ::Type{T}=M) where {M, T}

Constructs and returns an UnsteadyAeroInputs (ComponentVector) matching the model dimensions.
"""
function InputCache(model::UnsteadyAeroModel2D{M}, ::Type{T}=M) where {M, T}
    nCtrl = length(model.controlSurfaces)
    nVert = model.sizes.totalVertices
    
    return ComponentVector{T}(
        vb = zeros(T, 3),
        ab = zeros(T, 3),
        ωb = zeros(T, 3),
        dωb = zeros(T, 3),
        δc = zeros(T, nCtrl),
        dδc = zeros(T, nCtrl),
        ddδc = zeros(T, nCtrl),
        rs = zeros(T, 3 * nVert),
        vs = zeros(T, 3 * nVert),
        as = zeros(T, 3 * nVert),
        ρ = T[1.225]
    )
end

"""
    OutputCache(model::UnsteadyAeroModel2D{M}, ::Type{T}=M) where {M, T}

Constructs and returns an UnsteadyAeroOutputs (ComponentVector) matching the model dimensions.
"""
function OutputCache(model::UnsteadyAeroModel2D{M}, ::Type{T}=M) where {M, T}
    nmp = length(model.monitorPoints)
    return ComponentVector{T}(
        Fb = zeros(T, 3),
        Mb = zeros(T, 3),
        coeffsBody = zeros(T, 6),
        coeffsStab = zeros(T, 6),
        Fmp = zeros(T, 6 * nmp)
    )
end

# ==============================================================================
# 2. UnsteadyAeroCache
# ==============================================================================

"""
    UnsteadyAeroCache{T}

Standard internal cache struct for unsteady simulations.
"""
struct UnsteadyAeroCache{T}
    rVertex::Vector{Point3{T}}
    vVertex::Vector{Point3{T}}
    aVertex::Vector{Point3{T}}
    b::Vector{T}
    db::Vector{T}
    Γb::Vector{T}
    dΓb::Vector{T}
    Γw::Vector{T}
    Γs::Vector{T}
    Fa::Vector{Point3{T}}
    Ra::Vector{Point3{T}}
end

# --- Constructor ---
"""
    StateCache(model::UnsteadyAeroModel2D{M}, ::Type{T}=M) where {M, T}

Constructs and returns an UnsteadyAeroCache object matching the model dimensions.
"""
function StateCache(model::UnsteadyAeroModel2D{M}, ::Type{T}=M) where {M, T}
    nVert = model.sizes.totalVertices
    nPan = model.sizes.totalPanels
    nWake = model.wakeSizes.totalPanels
    nSSeg = model.sizes.totalSpanSegments
    nCSeg = model.sizes.totalChordSegments
    
    TArray = Point3{T}
    
    rVertex = zeros(TArray, nVert)
    vVertex = zeros(TArray, nVert)
    aVertex = zeros(TArray, nVert)
    b       = zeros(T, nPan)
    db      = zeros(T, nPan)
    Γb      = zeros(T, nPan)
    dΓb     = zeros(T, nPan)
    Γw      = zeros(T, nWake)
    Γs      = zeros(T, nSSeg + nCSeg)
    Fa      = zeros(TArray, nSSeg + nCSeg + nPan)
    Ra      = zeros(TArray, nSSeg + nCSeg + nPan)
    
    Ra[1:nSSeg + nCSeg] .= model.segmentProps.mid
    Ra[nSSeg + nCSeg + 1:end] .= model.panelProperties.rMid
    
    return UnsteadyAeroCache(
        rVertex, vVertex, aVertex, b, db,
        Γb, dΓb, Γw, Γs, Fa, Ra
    )
end
