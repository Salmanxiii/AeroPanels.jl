# Γw is the wake circulation solved from transport equation
# Γw0 is the leading edge circulation at the start of the wake which is equal to the body TE circulation i.e. Kutta condition
# Γb is the body circulation solved from AIC
#    Γw = K8 \ (-(K9*b)) (for steady case)
#    Γb = L3*Γw - L4*b
#    Γw0 = L7*Γw - L8*b

"""
$(SIGNATURES)

Calculate the unsteady aerodynamic influence coefficients and state-space matrices.
Based on Binder (2017) and AIAA (2018) publications.
"""
function UnsteadyWakeInfluence(rCollocation::Vector{Point3{T}}, normal, ringMesh, wakeMesh, bodySizes, wakeSizes, symmXZ, Δxw) where T
    AICwake = Influence(rCollocation, normal, wakeMesh, symmXZ);
    wakeLEIndices = LEPanelIndex(wakeSizes)
    wakeNonKuttaIndices = NonKuttaPanelIndex(wakeSizes)

    # Nonpermeability Condition
    # K1*Γb + K2*Γw0 + K3*Γw1 - b = 0
    K1 = Influence(rCollocation, normal, ringMesh, symmXZ);
    K2 = AICwake[:, wakeLEIndices]
    # K3 = @view AICwake[:, 1:end .∉ [wakeLEIndices]]
    K3 = @view AICwake[:, wakeNonKuttaIndices]

    bodyTEIndices = TEPanelIndex(bodySizes)
    totalKuttaPanels = length(wakeLEIndices)
    
    # Kutta Condition
    # K4*Γb + K5*Γw0 = 0
    K4 = SelectionOperator(bodyTEIndices, bodySizes.totalPanels,
        1:totalKuttaPanels, totalKuttaPanels)
    K5 = -I

    wSizes2 = Sizes([(nc-1, ns) for (s, nc, ns) in wakeSizes])
    totalWakePanels = wSizes2.totalPanels

    # Transport of Wake
    # dΓw1 = K6*Γw1 + K7*Γw0
    fromPanels = [PanelIndex(s,i,j, wSizes2) for (s, nc, ns) in wSizes2 for i in 1:nc-1 for j in 1:ns]
    toPanels = [PanelIndex(s,i,j, wSizes2) for (s, nc, ns) in wSizes2 for i in 2:nc for j in 1:ns]
    K6 = SelectionOperator(fromPanels, totalWakePanels,
    toPanels, totalWakePanels) - I
    K6 *= 1/Δxw 
    # Adds Kutta Panel Circulation to its trailing panels in wake
    K7 = SelectionOperator(1:totalKuttaPanels, totalKuttaPanels,
        LEPanelIndex(wSizes2), totalWakePanels)  
    K7 *= 1/Δxw

#  Circulation at wake LE panels is given in paper as :
#           Γw0 = L7*Γw - L8*b      (1)
#  But from kutta condition it is also equal to body TE circulation i.e. Kutta condition and circulation on body is given as:
#           Γb  = L3*Γw - L4*b      (2)
#  So equation (1) can also be written as:
#           Γw0  = Γb[TEindices] = L3[TEindices,:]*Γw - L4[TEindices,:]*b
# which means that:
#           L7 = L3[TEindices,:], L8 = L4[TEindices,:]
#  Thus L7 and L8 can be computed from L3 and L4, which is computationally more efficient.
#  Paper approach is left here for reference.
#  L8 = inv(K5 - K4K1⁻¹*K2) * K4K1⁻¹, L7 = L8 * K3
    
    # Γb = L3*Γw1 - L4*b
    L4 = inv(K2*K5*K4 - K1)
    L3 = L4 * K3
    
    # Γw0 = Γb[TEindices] = L7*Γw1 - L8*b
    L7 = L3[bodyTEIndices, :] #L7 = K4 * L3
    L8 = L4[bodyTEIndices, :] #L8 = K4 * L4
    
    # Γw1 = K8*Γw1 + K9*b
    # K9 should be sparse because K7 is highly sparse
    K9 = sparse(K7 * L8)
    K8 = sparse(K6 + K9*K3)
    K9 *= -1

    # dΓb = L5*Γw1 + L6*b - L4*db
    L5 = L3*K8
    L6 = L3*K9 # typo in paper 

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
struct UnsteadyAeroModel2D{T} <: AeroModel
    mesh::Mesh{3, T, QuadFace{Int64}, (:position,), Tuple{Vector{Point{3, T}}}, Vector{QuadFace{Int64}}}
    ringMesh::Mesh{3, T, QuadFace{Int64}, (:position,), Tuple{Vector{Point{3, T}}}, Vector{QuadFace{Int64}}}
    wakeMesh::Mesh{3, T, QuadFace{Int64}, (:position,), Tuple{Vector{Point{3, T}}}, Vector{QuadFace{Int64}}}
    sizes::Sizes
    wakeSizes::Sizes
    panelProperties::PanelProperties{T}
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
end

function UnsteadyAeroModel2D(surfaces::Vector{AeroSurface2D{T}}, props::AeroModelProperties{T}, V::T;
     nWake=80, wakeLength=20.,
     controlSurfaces=ControlSurface{T}[],
     monitorPoints=MonitorPoint{T}[]) where T

    mesh, sizes = CreateAeroMesh(surfaces)
    ringMesh = RingMesh(mesh, sizes)
    wakeMesh, wakeSizes = FlatWakeMesh(ringMesh, sizes, props; nWake=nWake, wakeLength=wakeLength)
    panelProperties = PanelProperties(sizes.totalPanels, mesh, FlowAxis(props))
    Δxw = wakeLength*props.c / nWake

    (K8, K9), (L3, L4, L5, L6, Γw0Indices, Γw1Indices, ΓbTEIndices) = UnsteadyWakeInfluence(panelProperties.rCollocation,
    panelProperties.normal, ringMesh, wakeMesh, sizes, wakeSizes, props.symmXZ, Δxw);
    segmentProps = ProcessSegments(ringMesh, sizes, wakeMesh, wakeSizes, props.symmXZ)
    
    return UnsteadyAeroModel2D(mesh, ringMesh, wakeMesh, sizes, wakeSizes, panelProperties,
    props, segmentProps, K8, K9, L3, L4, L5, L6, Γw0Indices, Γw1Indices, ΓbTEIndices, controlSurfaces, monitorPoints)
end

"""
    (model::UnsteadyAeroModel2D)(dΓw1, Γw1, p, t=0.0)

Functor for OrdinaryDiffEq. Solves for the wake transport states derivative.
`p` can be `b` (normalwash) or a tuple `(b, V)`.
"""
function (model::UnsteadyAeroModel2D)(dΓw1, Γw1, p, t=0.0)
    if p isa Tuple
        b, V = p
    else
        b, V = p, 1.0
    end
    mul!(dΓw1, model.K8, Γw1)
    mul!(dΓw1, model.K9, b, 1.0, 1.0)
    dΓw1 .*= V
    return nothing
end

"""
$(SIGNATURES)

Return the number of states in the unsteady state-space model.
"""
NumberOfStates(model::UnsteadyAeroModel2D) = size(model.K8, 1)
