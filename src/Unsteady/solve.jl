# --- Verification & Calculations Functions (Unmodified) ---

"""
    AddUnsteadyKinematics!(aVertex, rVertex, ab, dωb, ddδc, as, model)

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
    CalculateDNormalwash!(db, rVertex, vVertex, aVertex, model)

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
    SolveUnsteadyCirculation!(dΓw1, Γw1, dΓb, Γb, Γw, Γs, b, db, V, model::UnsteadyAeroModel2D)

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
    Γw[model.Γw1Indices] .= Γw1
    Γw[model.Γw0Indices] .= @view Γb[model.ΓbTEIndices]

    SegmentCirculation!(Γs, Γb, model.segmentProps)
    return nothing
end

"""
    CalculateUnsteadyAerodynamicForce!(Fa, dΓb, Γb, Γw, Γs, vVertex, model::UnsteadyAeroModel2D, ρ)

Calculate aerodynamic forces for the unsteady model.
"""
function CalculateUnsteadyAerodynamicForce!(Fa, dΓb, Γb, Γw, Γs, vVertex, model::UnsteadyAeroModel2D, ρ)
    nss = model.sizes.totalSpanSegments
    ncs = model.sizes.totalChordSegments
    nts = nss + ncs
    
    FaSteady = @view Fa[1:nts]
    CalculateAerodynamicForce!(FaSteady, Γb, Γw, Γs, vVertex, model, ρ)

    # Unsteady contribution: ρ * dΓ/dt * Area * Normal
    Fa_uns = @view Fa[nts+1:end]
    Fa_uns .= ρ .* dΓb .* model.panelProperties.area .* model.panelProperties.normal
    return nothing
end

# --- Split Solver Setup Functions ---

"""
    PrepareCacheForStates!(cache::UnsteadyAeroCache, inputs, model::UnsteadyAeroModel2D)

Calculates dynamic kinematics and panel normal wash in the cache.
"""
function PrepareCacheForStates!(cache::UnsteadyAeroCache, inputs, model::UnsteadyAeroModel2D)
    AddSteadyKinematics!(cache.rVertex, cache.vVertex, SVector{3}(inputs.vb), SVector{3}(inputs.ωb), 
                         inputs.δc, inputs.dδc, reinterpret(Point3{eltype(inputs)}, inputs.rs), reinterpret(Point3{eltype(inputs)}, inputs.vs), model)
    CalculateNormalwash!(cache.b, cache.rVertex, cache.vVertex, model)
    return nothing
end

"""
    PrepareCacheForOutputs!(cache::UnsteadyAeroCache, inputs, Γw1, dΓw1, model::UnsteadyAeroModel2D)

Calculates the unsteady wash derivative, wake circulations, and in-place aerodynamic forces.
"""
function PrepareCacheForOutputs!(cache::UnsteadyAeroCache, inputs, Γw1, dΓw1, model::UnsteadyAeroModel2D)
    AddUnsteadyKinematics!(cache.aVertex, cache.rVertex, SVector{3}(inputs.ab), SVector{3}(inputs.dωb), 
                           inputs.ddδc, reinterpret(Point3{eltype(inputs)}, inputs.as), model)
    CalculateDNormalwash!(cache.db, cache.rVertex, cache.vVertex, cache.aVertex, model)
    
    # Pre-allocated arrays inside cache receive the circulation derivative values
    V = norm(inputs.vb)
    SolveUnsteadyCirculation!(dΓw1, Γw1, cache.dΓb, cache.Γb, cache.Γw, cache.Γs, cache.b, cache.db, V, model)
    CalculateUnsteadyAerodynamicForce!(cache.Fa, cache.dΓb, cache.Γb, cache.Γw, cache.Γs, cache.vVertex, model, inputs.ρ[1])
    return nothing
end


# --- Time Integration Functor ---

function (model::UnsteadyAeroModel2D)(dΓw1, Γw1, p, t)
    input_func!, inputs, cache, _... = p

    # 1. Update inputs dynamically at intermediate solver stage time t
    input_func!(inputs, t)

    # 2. Recompute kinematics & normalwash to get new states coefficients
    PrepareCacheForStates!(cache, inputs, model)

    # 3. Calculate derivative
    mul!(dΓw1, model.K8, Γw1)
    mul!(dΓw1, model.K9, cache.b, 1.0, 1.0)
    V = norm(inputs.vb)
    dΓw1 .*= V
    return nothing
end

"""
    Outputs!(out, cache, model, vb, ρ)

Evaluates forces, moments, and stability/body coefficients in-place into `out`.
Uses standard coordinate conversion functions defined in `Misc.jl`.
"""
function Outputs!(out, cache, model, vb, ρ)
    Fg, Mg = GetTotalForces(cache, model)
    mProps = model.modelProperties
    
    # 1. Body axis forces and moments
    out.Fb .= GeometryToBodyAxis(Fg, mProps)
    out.Mb .= GeometryToBodyAxis(Mg, mProps)
    
    # Aerodynamic angles (angle of attack α)
    α, _, _ = AerodynamicAngles(vb)
    
    # 2. Body axis coefficients (unrolled direct assignment)
    CFbody, CMbody = ConvertToCoefficients(out.Fb, out.Mb, vb, ρ, mProps)
    out.coeffsBody .= (CFbody..., CMbody...)
    
    # 3. Stability axis coefficients (unrolled direct assignment)
    Fstab = GeometryToStabilityAxis(Fg, α, mProps)
    Mstab = GeometryToStabilityAxis(Mg, α, mProps)
    CFstab, CMstab = ConvertToCoefficients(Fstab, Mstab, vb, ρ, mProps)
    out.coeffsStab .= (CFstab..., CMstab...)
    
    # 4. Monitor point loads
    n_mp = length(model.monitorPoints)
    if n_mp > 0
        MonitorPointLoads!(out.Fmp, cache, model)
    end
    return nothing
end

# --- SavingCallback Save Function ---

# AeroSaveFunc(u, t, integrator)
#
# Callback function for SavingCallback. Updates inputs, prepares states, resolves
# circulation and forces, and returns a copied UnsteadyAeroOutputs structure.
function AeroSaveFunc(u, t, integrator)
    input_func!, inputs, cache, model, dΓw1_temp = integrator.p
    
    # 1. Update inputs
    input_func!(inputs, t)
    
    # 2. Prepare states (kinematics, b, V)
    PrepareCacheForStates!(cache, inputs, model)
    
    # 3. Evaluate state derivative in-place
    integrator.f(dΓw1_temp, u, integrator.p, t)
    
    # 4. Prepare outputs (unsteady wash, circulation derivatives, forces)
    PrepareCacheForOutputs!(cache, inputs, u, dΓw1_temp, model)
    
    # 5. Build outputs struct
    out = OutputCache(model, eltype(u))
    Outputs!(out, cache, model, SVector{3}(inputs.vb), inputs.ρ[1])
    return copy(out)
end

# --- User-Facing solver API ---

"""
    AeroSolve(model::UnsteadyAeroModel2D{T}, tspan, input_func!; solver, steady=false) where T

Solve an unsteady time-domain simulation using SciML SavingCallback.
"""
function AeroSolve(model::UnsteadyAeroModel2D{T}, tspan, input_func!; solver, steady=false) where T
    # 1. Allocate caches using the new constructor API
    inputs = InputCache(model, T)
    cache = StateCache(model, T)
    
    nWake1 = size(model.K8, 1)
    dΓw1_temp = zeros(T, nWake1)
    
    # ODE parameters tuple (input_func!, inputs, cache, model, dΓw1_temp)
    p = (input_func!, inputs, cache, model, dΓw1_temp)

    # 2. Handle optional steady initialization
    Γw1_0 = zeros(T, nWake1)
    if steady
        input_func!(inputs, 0.0)
        PrepareCacheForStates!(cache, inputs, model)
        Γw1_0 .= -model.K8 \ (model.K9 * cache.b)
    end
    
    # 3. Setup the ODE problem
    prob = ODEProblem(model, Γw1_0, tspan, p)
    
    # 4. Setup SavingCallback
    outputs = SavedValues(Float64, typeof(OutputCache(model, T)))
    scb = SavingCallback(AeroSaveFunc, outputs)
    
    # 5. Solve using high-level SciML solver with callback
    sol = solve(prob, solver; callback=scb)
    
    return (sol, outputs)
end