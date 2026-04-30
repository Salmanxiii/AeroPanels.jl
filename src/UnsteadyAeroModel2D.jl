# Γw is the wake circulation solved from transport equation
# Γw0 is the leading edge circulation at the start of the wake which is equal to the body TE circulation i.e. Kutta condition
# Γb is the body circulation solved from AIC
#    Γw = K8 \ (-(K9*b)) (for steady case)
#    Γb = L3*Γw - L4*b
#    Γw0 = L7*Γw - L8*b

function UnsteadyWakeInfluence(rCollocation::Vector{Point3{T}}, normal, ringMesh, wakeMesh, bodySizes, wakeSizes, symmXZ, Δt) where T
    AICwake = Influence(rCollocation, normal, wakeMesh, symmXZ);
    wakeLEIndices = LEPanelIndex(wakeSizes)
    wakeNonKuttaIndices = NonKuttaPanelIndex(wakeSizes)

    K1 = Influence(rCollocation, normal, ringMesh, symmXZ);
    K2 = AICwake[:, wakeLEIndices]
    # K3 = @view AICwake[:, 1:end .∉ [wakeLEIndices]]
    K3 = @view AICwake[:, wakeNonKuttaIndices]

    bodyTEIndices = TEPanelIndex(bodySizes)
    totalKuttaPanels = length(wakeLEIndices)
    # Kutta Convection
    K4 = SelectionOperator(bodyTEIndices, bodySizes.totalPanels,
        1:totalKuttaPanels, totalKuttaPanels)
    K5 = -I

    wSizes2 = Sizes([(nc-1, ns) for (s, nc, ns) in wakeSizes])
    totalWakePanels = wSizes2.totalPanels
    # Transport of Wake
    fromPanels = [PanelIndex(s,i,j, wSizes2) for (s, nc, ns) in wSizes2 for i in 1:nc-1 for j in 1:ns]
    toPanels = [PanelIndex(s,i,j, wSizes2) for (s, nc, ns) in wSizes2 for i in 2:nc for j in 1:ns]
    K6 = SelectionOperator(fromPanels, totalWakePanels,
    toPanels, totalWakePanels) - I
    # Adds Kutta Panel Circulation to its trailing panels in wake
    K7 = SelectionOperator(1:totalKuttaPanels, totalKuttaPanels,
        LEPanelIndex(wSizes2), totalWakePanels)
    K6 *= 1/Δt # K6 *= V / Δxw
    K7 *= 1/Δt
#  Circulation at wake LE panels is given in paper as :
#           Γw0 = L7*Γw - L8*b      (1)
#  But from kutta condition it is also equal to body TE panels and circulation on body is given as:
#           Γb  = L3*Γw - L4*b      (2)
#  So equation (1) can also be written as:
#           Γw0  = Γb[TEindices] = L3[TEindices,:]*Γw - L4[TEindices,:]*b
# which means that:
#           L7 = L3[TEindices,:], L8 = L4[TEindices,:]
#  Thus L7 and L8 can be computed from L3 and L4, which is computationally more efficient.
#  Paper approach is left here for reference.
#  L8 = inv(K5 - K4K1⁻¹*K2) * K4K1⁻¹, L7 = L8 * K3

    L4 = -inv(K1-K2*K5*K4)
    L8 = L4[bodyTEIndices, :]
    # K9 should be sparse because K7 is highly sparse
    #K4K1⁻¹ = K4 / K1
    #K_temp1 = inv(K5 - K4K1⁻¹*K2)
    #L8 = K_temp1 * K4K1⁻¹
    K9 = K7 * L8
    K8 = K6 + K9*K3
    K9 *= -1

    #L4 = -inv(K1-K2*K5*K4)
    L3 = L4*K3
    L5 = L3*K8
    L6 = L3*K9 # typo in paper 
    #L7 = L8 * K3
    # L8 calculated above before K9

    L7 = L3[bodyTEIndices, :]
    L9, L10 = FullWakeFromTransportWakeOperator(bodySizes, wakeSizes, L7, L8)

    return (sparse(K8), sparse(K9)), (L3, L4, L5, L6, L7, L8, L9, L10)
end

struct UnsteadyAeroModel2D{T} <: AeroModel
    mesh::GeometryBasics.Mesh{3, T}
    ringMesh::GeometryBasics.Mesh{3, T}
    wakeMesh::GeometryBasics.Mesh{3, T}
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
    L7::Matrix{T}
    L8::Matrix{T}
    L9::SparseMatrixCSC{T, Int}
    L10::SparseMatrixCSC{T, Int}
end

function UnsteadyAeroModel2D(surfaces::Vector{AeroSurface2D{T}}, props::AeroModelProperties{T}, V::T; nWake=80, wakeLength=20.) where T
    mesh, sizes = CreateAeroMesh(surfaces)
    ringMesh = RingMesh(mesh, sizes)
    wakeMesh, wakeSizes = FlatWakeMesh(ringMesh, sizes, props; nWake=nWake, wakeLength=wakeLength)
    panelProperties = PanelProperties(sizes.totalPanels, mesh, FlowAxis(props))
    Δxw = wakeLength*props.c / nWake
    Δt = Δxw / V

    (K8, K9), (L3, L4, L5, L6, L7, L8, L9, L10) = UnsteadyWakeInfluence(panelProperties.rCollocation,
    panelProperties.normal, ringMesh, wakeMesh, sizes, wakeSizes, props.symmXZ, Δt);
    segmentProps = ProcessSegments(ringMesh, sizes, wakeMesh, wakeSizes, props.symmXZ)
    return UnsteadyAeroModel2D(mesh, ringMesh, wakeMesh, sizes, wakeSizes, panelProperties,
    props, segmentProps, K8, K9, L3, L4, L5, L6, L7, L8, L9, L10)
end

function FullWakeFromTransportWakeOperator(bodySizes, wakeSizes, L7, L8)
    # Γw0(LE) = L7*Γw(transport) - L8*b
    # Γwake(full wake) = L9 * Γw + L10 * b

    wakeLEIndices = LEPanelIndex(wakeSizes)
    wakeNonKuttaIndices = NonKuttaPanelIndex(wakeSizes)
    wakeTotal = wakeSizes.totalPanels

    # Construct sparse operator
    m1, m2  = length(wakeLEIndices), length(wakeNonKuttaIndices)
    n = m1*m2 + length(wakeNonKuttaIndices)
    vals = ones(n)
    rows, cols = ones(Int, n), ones(Int, n)
    n=1
    for i in 1:m1
        for j in 1:m2
            rows[n] = wakeLEIndices[i]
            cols[n] = j
            n+=1
        end
    end
    m=1
    for i in 1:length(wakeNonKuttaIndices)
        rows[n] = wakeNonKuttaIndices[i]
        cols[n] = m
        n+=1
        m+=1
    end
    L9 =  sparse(rows, cols, vals, wakeTotal, m2)
    L9[wakeLEIndices, :] .= L7
    
    bodyTotal = bodySizes.totalPanels
    rows = [wakeLEIndices[j] for i in 1:bodyTotal for j in 1:m1]
    cols = [i for i in 1:bodyTotal for j in 1:m1]
    vals = vec(L8)
    vals *= -1
    L10 = sparse(rows, cols, vals, wakeTotal, bodyTotal) 

    return L9, L10
end

function SteadyForce(Γw, b, vb, ρ, model)
    Γb = model.L3*Γw - model.L4*b
    Γwake = GetFullWakeVector(Γw, b, model)
    Γs = SegmentCirculation(Γb, model.segmentProps)
    Fa = AerodynamicForce(Γb, Γwake, Γs, vb, model, ρ=ρ)
    return Fa, Γb
end

function SolveSteadyCirculation(b, model::UnsteadyAeroModel2D)
    Γw = - (model.K8 \ (model.K9*b))
    return Γw
end

function GetFullWakeVector(Γw::Vector, b::Vector, model::UnsteadyAeroModel2D)
    Γwake = model.L9 * Γw + model.L10 * b
    return Γwake
end

function AeroSolve(vb, aeroModel::UnsteadyAeroModel2D, ρ = 1.225)
    b = NormalWash(vb, aeroModel)
    Γw = SolveSteadyCirculation(b, aeroModel)
    Fa, _ = SteadyForce(Γw, b, vb, ρ, aeroModel)
    sol = SteadySolution(Fa, vb, aeroModel, ρ)
    return sol
end

function SolveCirculation(Γw::Vector, vb, model::UnsteadyAeroModel2D)
    b = NormalWash(vb, model)
    dΓw = SolveCirculation(Γw, b, model)
    return dΓw
end

function SolveCirculation(Γw::Vector, b::Vector, model::UnsteadyAeroModel2D)
    dΓw = model.K8*Γw + model.K9*b
    return dΓw
end

"""
Returns unsteady aerodynamic forces at panel center
"""
function UnsteadyPanelForces(Γw, b, ρ, model, db)
    dΓb = model.L5*Γw .+ model.L6*b .- model.L4*db
    Fu = ρ .* dΓb .* model.panelProperties.area .* model.panelProperties.normal
    return Fu, dΓb
end

function UnsteadyPanelForces(Γw, b, ρ, model)
    dΓb = model.L5*Γw .+ model.L6*b
    Fu = ρ .* dΓb .* model.panelProperties.area .* model.panelProperties.normal
    return Fu, dΓb
end

"""
Returns total aerodynamic forces including steady and unsteady contributions
"""
function SolveForces(Γw, vb, dvb, model::UnsteadyAeroModel2D, ρ = 1.225)
    b = NormalWash(vb, model)
    db = NormalWash(dvb, model)
    Fu, _ = UnsteadyPanelForces(Γw, b, ρ, model, db)
    Fa, _ = SteadyForce(Γw, b, vb, ρ, model)
    return Fa, Fu
end

NumberOfStates(model::UnsteadyAeroModel2D) = size(model.K8, 1)