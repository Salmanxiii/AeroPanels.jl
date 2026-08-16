# Γw is the wake circulation solved from transport equation
# Γw0 is the leading edge circulation at the start of the wake which is equal to the body TE circulation i.e. Kutta condition
# Γb is the body circulation solved from AIC
#    Γw = K8 \ (-(K9*b)) (for steady case)
#    Γb = L3*Γw - L4*b
#    Γw0 = L7*Γw - L8*b

"""
$(SIGNATURES)

Calculate the unsteady aerodynamic influence coefficients and state-space matrices.
Based on Binder (2017) and Werter (2018) publications.
"""
function UnsteadyWakeInfluence(rCollocation::Vector{Point3{T}}, normal, ringMesh, wakeMesh, bodySizes, wakeSizes, symmXZ, wakeSpacing::Vector{T}) where T
    AICwake = Influence(rCollocation, normal, wakeMesh, symmXZ);
    wakeLEIndices = LEPanelIndex(wakeSizes)
    wakeNonKuttaIndices = NonKuttaPanelIndex(wakeSizes)

    K1 = Influence(rCollocation, normal, ringMesh, symmXZ);
    K2 = AICwake[:, wakeLEIndices]
    K3 = AICwake[:, wakeNonKuttaIndices]

    bodyTEIndices = TEPanelIndex(bodySizes)
    totalKuttaPanels = length(wakeLEIndices)
    
    K4 = SelectionOperator(bodyTEIndices, bodySizes.totalPanels,
        1:totalKuttaPanels, totalKuttaPanels)
    K5 = -I

    wSizes2 = Sizes([(nc-1, ns) for (s, nc, ns) in wakeSizes])
    totalWakePanels = wSizes2.totalPanels

    fromPanels = [PanelIndex(s,i,j, wSizes2) for (s, nc, ns) in wSizes2 for i in 1:nc-1 for j in 1:ns]
    toPanels = [PanelIndex(s,i,j, wSizes2) for (s, nc, ns) in wSizes2 for i in 2:nc for j in 1:ns]
    K6 = SelectionOperator(fromPanels, totalWakePanels,
    toPanels, totalWakePanels) - I
    
    K7 = SelectionOperator(collect(1:totalKuttaPanels), totalKuttaPanels,
        LEPanelIndex(wSizes2), totalWakePanels)  
        
    inv_dx = zeros(T, totalWakePanels)
    for (s, nc, ns) in wSizes2
        for i in 1:nc
            for j in 1:ns
                index = PanelIndex(s, i, j, wSizes2)
                inv_dx[index] = 1.0 / wakeSpacing[i + 1]
            end
        end
    end
    
    D_inv_dx = sparse(Diagonal(inv_dx))
    K6 = D_inv_dx * K6
    K7 = D_inv_dx * K7
    
    L4 = inv(K2*K5*K4 - K1)
    L3 = L4 * K3
    
    L7 = L3[bodyTEIndices, :]
    L8 = L4[bodyTEIndices, :]
    
    K9 = sparse(K7 * L8)
    K8 = sparse(K6 + K9*K3)
    K9 *= -1

    L5 = L3*K8
    L6 = L3*K9

    Γw0Indices = LEPanelIndex(wakeSizes)
    Γw1Indices = NonKuttaPanelIndex(wakeSizes)
    ΓbTEIndices = bodyTEIndices

    return (K8, K9), (L3, L4, L5, L6, Γw0Indices, Γw1Indices, ΓbTEIndices)
end

"""
$(TYPEDEF)

An unsteady 2D aerodynamic model based on the Continuous-Time Unsteady Vortex Lattice Method.

$(TYPEDFIELDS)
"""
@fmu_model struct UnsteadyAeroModel2D{T} <: AbstractFMUModel{T}
    boundMesh::ThinAeroMesh{T}
    wakeMesh::FixedWakeMesh{T}
    modelProperties::AeroModelProperties{T}
    segmentProps::SegmentProperties{T}
    K8::SparseMatrixCSC{T, Int}
    K9::SparseMatrixCSC{T, Int}
    L3::Matrix{T}
    L4::Matrix{T}
    L5::Matrix{T}
    L6::Matrix{T}
    Γw0Indices::Vector{Int}
    Γw1Indices::Vector{Int}
    ΓbTEIndices::Vector{Int}
    controlSurfaces::Vector{ControlSurface{T}}
    monitorPoints::Vector{MonitorPoint{T}}
    
    [States]
    circ_w1 = (m) -> size(m.K8, 1)

    [Derivatives]
    dcirc_w1 = (m) -> size(m.K8, 1)

    [Inputs]
    vb     = 3
    wb     = 3
    ab     = 3
    dwb    = 3
    u_cs   = (m) -> length(m.controlSurfaces)
    du_cs  = (m) -> length(m.controlSurfaces)
    ddu_cs = (m) -> length(m.controlSurfaces)
    rs     = (m) -> 3 * m.boundMesh.sizes.totalVertices
    vs     = (m) -> 3 * m.boundMesh.sizes.totalVertices
    as     = (m) -> 3 * m.boundMesh.sizes.totalVertices
    rho    = 1

    [Outputs]
    Fb          = 3
    Mb          = 3
    Fb_unsteady = 3
    Cb          = 6
    Cs          = 6
    Fmp         = (m) -> 6 * length(m.monitorPoints)
    Fa          = (m) -> 3 * (TotalSegments(m.boundMesh) + TotalPanels(m.boundMesh))
end

function UnsteadyAeroModel2D(boundMesh::ThinAeroMesh{T}, wakeMesh::FixedWakeMesh{T}, props::AeroModelProperties{T}, V::T;
     controlSurfaces=ControlSurface{T}[],
     monitorPoints=MonitorPoint{T}[]) where T

    wakeSpacing = wakeMesh.wakeSpacing

    rCollocation = GetCollocationPoints(boundMesh)
    normals = GetNormalVectors(boundMesh)
    
    (K8, K9), (L3, L4, L5, L6, Γw0Indices, Γw1Indices, ΓbTEIndices) = UnsteadyWakeInfluence(rCollocation,
    normals, boundMesh.ringMesh, wakeMesh.mesh, boundMesh.sizes, wakeMesh.wakeSizes, props.symmXZ, wakeSpacing);
    segmentProps = ProcessSegments(boundMesh.ringMesh, boundMesh.sizes, wakeMesh.mesh, wakeMesh.wakeSizes, props.symmXZ)

    return UnsteadyAeroModel2D(boundMesh, wakeMesh, props, segmentProps, K8, K9, L3, L4, L5, L6, Γw0Indices, Γw1Indices, ΓbTEIndices, controlSurfaces, monitorPoints)
end



"""
    UnsteadyAeroModel2DCache{T}

Monolithic internal cache struct for zero-allocation unsteady evaluations.
"""
struct UnsteadyAeroModel2DCache{T}
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

function AeroPanels.AllocateFMUCaches(model::UnsteadyAeroModel2D{M}, ::Type{T}) where {M, T}
    nVert = model.boundMesh.sizes.totalVertices
    nPan = model.boundMesh.sizes.totalPanels
    nWake = model.wakeMesh.wakeSizes.totalPanels
    nSSeg = model.boundMesh.sizes.totalSpanSegments
    nCSeg = model.boundMesh.sizes.totalChordSegments
    
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
    Ra[nSSeg + nCSeg + 1:end] .= GetMidpoints(model.boundMesh)
    
    return UnsteadyAeroModel2DCache(
        rVertex, vVertex, aVertex, b, db,
        Γb, dΓb, Γw, Γs, Fa, Ra
    )
end
