"""
$(TYPEDEF)

A struct holding the cache arrays for unsteady aerodynamic simulations to minimize allocations.

$(TYPEDFIELDS)
"""
struct UnsteadyAeroCache{T}
    rVertex::Vector{Point3{T}}
    vVertex::Vector{Point3{T}}
    aVertex::Vector{Point3{T}}
    b::Vector{T}
    db::Vector{T}
    Γw1::Vector{T}
    dΓw1::Vector{T}
    Γb::Vector{T}
    dΓb::Vector{T}
    Γw::Vector{T}
    Γs::Vector{T}
    Fa::Vector{Point3{T}}
    Ra::Vector{Point3{T}}
end

"""
$(SIGNATURES)

Create and return a cache of arrays for the unsteady solver.
"""
function CreateCacheArrays(model::UnsteadyAeroModel2D{T}) where T
    return CreateCacheArrays(model, T)
end

function CreateCacheArrays(model::UnsteadyAeroModel2D, T)
    nVert = model.sizes.totalVertices
    nPan = model.sizes.totalPanels
    nWake1 = NumberOfStates(model)
    nWake = model.wakeSizes.totalPanels
    nSSeg = model.sizes.totalSpanSegments
    nCSeg = model.sizes.totalChordSegments

    TArray = Point3{T}
    rVertex = zeros(TArray, nVert)
    vVertex = zeros(TArray, nVert)
    aVertex = zeros(TArray, nVert)
    b = zeros(T, nPan)
    db = zeros(T, nPan)
    Γw1 = zeros(T, nWake1)
    dΓw1 = zeros(T, nWake1)
    Γb = zeros(T, nPan)
    dΓb = zeros(T, nPan)
    Γw = zeros(T, nWake)
    Γs = zeros(T, nSSeg + nCSeg)
    Fa = zeros(TArray, nSSeg + nCSeg + nPan)
    Ra = zeros(TArray, nSSeg + nCSeg + nPan)
    Ra[1:nSSeg + nCSeg] .= model.segmentProps.mid
    Ra[nSSeg + nCSeg + 1:end] .= model.panelProperties.rMid

    return UnsteadyAeroCache(rVertex, vVertex, aVertex, b, db, 
            Γw1, dΓw1, Γb, dΓb, Γw, Γs, Fa, Ra)
end

"""
$(SIGNATURES)

Calculate the unsteady kinematics (vertex accelerations).
"""
function AddUnsteadyKinematics!(aVertex, rVertex,
    ab, dωb, ddδc, as, model)

    fill!(aVertex, zero(eltype(aVertex)))
    CSAcceleration!(aVertex, ddδc, model)

    bodyAxis = model.modelProperties.bodyFixedCS
    CG = model.modelProperties.CG
    (ab_g, dωb_g) = (bodyAxis' * component for component in (ab, dωb))

    for i in eachindex(aVertex)
        aVertex[i] = aVertex[i] - ab_g + as[i] - cross(dωb_g, rVertex[i] - CG)
    end
    return nothing
end

"""
$(SIGNATURES)

Calculate the time derivative of the normal wash vector `db`.
"""
function CalculateDNormalwash!(db, rVertex, vVertex, aVertex, model)
    # db = - (a . n + v . dn)
    for (s, nc, ns) in model.sizes
        @batch for i in 1:nc
            for j in 1:ns
                p = PanelIndex(s, i, j, model.sizes)
                i1, i2, i3, i4 = PanelVertexIndices(s, i, j, model.sizes)
                r1, r2, r3, r4 = rVertex[i1], rVertex[i2], rVertex[i3], rVertex[i4]
                v1, v2, v3, v4 = vVertex[i1], vVertex[i2], vVertex[i3], vVertex[i4]
                a1, a2, a3, a4 = aVertex[i1], aVertex[i2], aVertex[i3], aVertex[i4]

                normalVec = cross(r3 - r1, r2 - r4)
                normalNorm = norm(normalVec)
                normalUnitVec = normalVec / normalNorm
                
                v = 0.125*(v1 + v2) + 0.375*(v3 + v4)
                a = 0.125*(a1 + a2) + 0.375*(a3 + a4)

                dNormalVec = cross(v3 - v1, r2 - r4) + cross(r3 - r1, v2 - v4)
                dNormalUnitVec = (dNormalVec - normalUnitVec * dot(normalUnitVec, dNormalVec)) / normalNorm
                db[p] = -(dot(a, normalUnitVec) + dot(v, dNormalUnitVec))
            end
        end
    end
    return nothing
end

"""
$(SIGNATURES)

Solve for the time derivatives and intermediate circulations of the unsteady system.
"""
function SolveUnsteadyCirculation!(dΓw1, Γw1, dΓb, Γb, Γw, Γs, b, db, V, model::UnsteadyAeroModel2D)
    # dΓw1 = V * (K8*Γw1 + K9*b)
    mul!(dΓw1, model.K8, Γw1)
    mul!(dΓw1, model.K9, b, 1.0, 1.0)
    dΓw1 .*= V

    # Γb = L3*Γw1 - L4*b
    mul!(Γb, model.L3, Γw1)
    mul!(Γb, model.L4, b, -1.0, 1.0)

    # dΓb = V * (L5*Γw1 + L6*b) - L4*db
    mul!(dΓb, model.L5, Γw1)
    mul!(dΓb, model.L6, b, 1.0, 1.0)
    dΓb .*= V
    mul!(dΓb, model.L4, db, -1.0, 1.0)

    # Reconstruct the full wake circulation vector Γw from transport states and Kutta condition
    # Γw[transport_indices] = Γw1
    # Γw[kutta_indices] = Γb[TE_indices]
    Γw[model.Γw1Indices] .= Γw1
    Γw[model.Γw0Indices] .= @view Γb[model.ΓbTEIndices]

    SegmentCirculation!(Γs, Γb, model.segmentProps)
    return nothing
end

"""
$(SIGNATURES)

Calculate aerodynamic forces for the unsteady model.
"""
function CalculateUnsteadyAerodynamicForce!(Fa, dΓb, Γb, Γw, Γs, vVertex, model::UnsteadyAeroModel2D, ρ)
    nss = model.sizes.totalSpanSegments
    ncs = model.sizes.totalChordSegments
    nts = nss + ncs
    n = model.sizes.totalPanels
    
    FaSteady = @view Fa[1:nts]
    CalculateAerodynamicForce!(FaSteady, Γb, Γw, Γs, vVertex, model, ρ)

    # Unsteady contribution: ρ * dΓ/dt * Area * Normal
    @batch for i in 1:n
        Fa[nts+i] = ρ * dΓb[i] * model.panelProperties.area[i] * model.panelProperties.normal[i]
    end
    return nothing
end

"""
$(SIGNATURES)

In-place version of the unsteady `AeroSolve`.
"""
function AeroSolve!(cache, vb, ab, ωb, dωb, δc, dδc, ddδc, rs, vs, as, Γw1, model::UnsteadyAeroModel2D, ρ=1.225)
    AddSteadyKinematics!(cache.rVertex, cache.vVertex, vb, ωb, δc, dδc, rs, vs, model)
    AddUnsteadyKinematics!(cache.aVertex, cache.rVertex, ab, dωb, ddδc, as, model)
    
    CalculateNormalwash!(cache.b, cache.rVertex, cache.vVertex, model)
    CalculateDNormalwash!(cache.db, cache.rVertex, cache.vVertex, cache.aVertex, model)
    
    V = norm(vb)
    SolveUnsteadyCirculation!(cache.dΓw1, Γw1, cache.dΓb, cache.Γb, cache.Γw, cache.Γs, cache.b, cache.db, V, model)
    CalculateUnsteadyAerodynamicForce!(cache.Fa, cache.dΓb, cache.Γb, cache.Γw, cache.Γs, cache.vVertex, model, ρ)
    return nothing
end

"""
    AeroSolve(model::UnsteadyAeroModel2D, Γw1_0, tspan, input_func!; solver)

Solve an unsteady time-domain simulation and auto-compute forces/coefficients.
"""
function AeroSolve(model::UnsteadyAeroModel2D, Γw1_0, tspan, input_func!; solver)
    T = eltype(Γw1_0)
    inputs = AeroInputs(model, T=T)
    cache = CreateCacheArrays(model, T)
    
    p = UnsteadySimParams(input_func!, inputs, cache)
    prob = ODEProblem(model, Γw1_0, tspan, p)
    
    sol = solve(prob, solver)
    
    n_steps = length(sol.t)
    CX = zeros(T, n_steps)
    CY = zeros(T, n_steps)
    CZ = zeros(T, n_steps)
    CD = zeros(T, n_steps)
    CL = zeros(T, n_steps)
    F_tot = Vector{SVector{3, T}}(undef, n_steps)
    M_tot = Vector{SVector{3, T}}(undef, n_steps)
    
    n_mp = length(model.monitorPoints)
    mp_loads = [Vector{SVector{6, T}}(undef, n_steps) for _ in 1:n_mp]
    Fmp = zeros(T, 6 * max(1, n_mp)) # Preallocate
    
    # Post-process for all time steps
    for (i, t) in enumerate(sol.t)
        input_func!(inputs, t)
        Γw1 = sol.u[i]
        
        AeroSolve!(cache, SVector(inputs.vb), SVector(inputs.ab), SVector(inputs.ωb), SVector(inputs.dωb), 
                   inputs.δc, inputs.dδc, inputs.ddδc, inputs.rs, inputs.vs, inputs.as, Γw1, model, inputs.ρ)
        
        F, M = GetTotalForces(cache, model)
        F_tot[i] = F
        M_tot[i] = M
        
        CFstab, CMstab = GetStabilityCoefficients(cache, SVector(inputs.vb), inputs.ρ, model)
        CX[i] = CFstab[1]
        CY[i] = CFstab[2]
        CZ[i] = CFstab[3]
        CD[i] = CFstab[1]
        CL[i] = -CFstab[3]
        
        if n_mp > 0
            MonitorPointLoads!(Fmp, cache, model)
            for j in 1:n_mp
                mp_loads[j][i] = SVector{6, T}(Fmp[6*(j-1)+1 : 6*(j-1)+6]...)
            end
        end
    end
    
    return (t=sol.t, CX=CX, CY=CY, CZ=CZ, CD=CD, CL=CL, F=F_tot, M=M_tot, mp_loads=mp_loads, sol=sol)
end


"""
$(TYPEDEF)

A struct holding the time-dependent inputs for unsteady aerodynamic simulations.
These inputs are meant to be updated in-place during time integration.
"""
mutable struct AeroInputs{T}
    vb::MVector{3, T}
    ab::MVector{3, T}
    ωb::MVector{3, T}
    dωb::MVector{3, T}
    δc::Vector{T}
    dδc::Vector{T}
    ddδc::Vector{T}
    rs::Vector{Point3{T}}
    vs::Vector{Point3{T}}
    as::Vector{Point3{T}}
    ρ::T
end

function AeroInputs(model::UnsteadyAeroModel2D; T=Float64)
    nVert = model.sizes.totalVertices
    nCtrl = length(model.controlSurfaces)
    return AeroInputs{T}(
        MVector{3, T}(0,0,0),
        MVector{3, T}(0,0,0),
        MVector{3, T}(0,0,0),
        MVector{3, T}(0,0,0),
        zeros(T, nCtrl),
        zeros(T, nCtrl),
        zeros(T, nCtrl),
        fill(zero(Point3{T}), nVert),
        fill(zero(Point3{T}), nVert),
        fill(zero(Point3{T}), nVert),
        1.225
    )
end

function AeroInputs(model::UnsteadyAeroModel2D{Float64})
    nVert = model.sizes.totalVertices
    nCtrl = length(model.controlSurfaces)
    return AeroInputs{Float64}(
        MVector{3, Float64}(0,0,0),
        MVector{3, Float64}(0,0,0),
        MVector{3, Float64}(0,0,0),
        MVector{3, Float64}(0,0,0),
        zeros(Float64, nCtrl),
        zeros(Float64, nCtrl),
        zeros(Float64, nCtrl),
        fill(zero(Point3{Float64}), nVert),
        fill(zero(Point3{Float64}), nVert),
        fill(zero(Point3{Float64}), nVert),
        1.225
    )
end


struct UnsteadySimParams{F, I, C}
    input_func!::F
    inputs::I
    cache::C
end

function (model::UnsteadyAeroModel2D)(dΓw1, Γw1, p::UnsteadySimParams, t)
    # 1. Update inputs for current t
    p.input_func!(p.inputs, t)

    # 2. Recompute steady kinematics to get new normalwash
    AddSteadyKinematics!(p.cache.rVertex, p.cache.vVertex, SVector(p.inputs.vb), SVector(p.inputs.ωb),
                         p.inputs.δc, p.inputs.dδc, p.inputs.rs, p.inputs.vs, model)
    CalculateNormalwash!(p.cache.b, p.cache.rVertex, p.cache.vVertex, model)

    V = norm(p.inputs.vb)

    # 3. Calculate derivative
    mul!(dΓw1, model.K8, Γw1)
    mul!(dΓw1, model.K9, p.cache.b, 1.0, 1.0)
    dΓw1 .*= V
    return nothing
end



