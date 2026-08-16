# src/Models/Dynamics/RigidBody6DOF.jl

"""
    AddedMass(aero::UnsteadyAeroModel2D{T}) where T

Compute the 6x6 coupled added mass and added inertia matrix at reference density rho = 1.0.
Returns a 6x6 matrix partition:
    [ M_added   C_added ]
    [ D_added   I_added ]
"""
function AddedMass(aero::UnsteadyAeroModel2D{T}) where T
    panels = aero.boundMesh
    Np = length(panels.areas)
    rRef = aero.modelProperties.rRef
    R = aero.modelProperties.bodyFixedCS
    
    # 1. W_full (normalwash kinematics) uses rCollocation
    W_full = Matrix{T}(undef, Np, 6)
    for i in 1:Np
        n_g = panels.normals[i]
        r_g = panels.rCollocation[i] - rRef
        n = R * n_g
        rxn = R * cross(r_g, n_g)
        
        W_full[i, 1] = n[1]; W_full[i, 2] = n[2]; W_full[i, 3] = n[3]
        W_full[i, 4] = rxn[1]; W_full[i, 5] = rxn[2]; W_full[i, 6] = rxn[3]
    end
    
    # 2. W_force (moment projection) uses rMid
    W_force = Matrix{T}(undef, Np, 6)
    for i in 1:Np
        n_g = panels.normals[i]
        r_g = panels.rMid[i] - rRef
        n = R * n_g
        rxn = R * cross(r_g, n_g)
        
        W_force[i, 1] = n[1]; W_force[i, 2] = n[2]; W_force[i, 3] = n[3]
        W_force[i, 4] = rxn[1]; W_force[i, 5] = rxn[2]; W_force[i, 6] = rxn[3]
    end
    
    # Pre-multiply panel areas into W_force
    W_area = panels.areas .* W_force
    
    # Compute 6x6 added mass matrix at reference density rho = 1.0
    return W_area' * (aero.L4 * W_full)
end

"""
    RigidBody6DOF{M, T}

Struct to represent a coupled 6DOF flight dynamics model with unsteady aerodynamics.
"""
@fmu_model struct RigidBody6DOF{M, T} <: AbstractFMUModel{T}
    aero_model::M
    mass::T
    inertia::SMatrix{3, 3, T, 9}
    inertia_inv::SMatrix{3, 3, T, 9}
    massAero::SMatrix{6, 6, T, 36}

    [States]
    pos      = (size=3, output=true)
    euler    = (size=3, output=true)
    vb       = (size=3, output=true)
    wb       = (size=3, output=true)
    circ_w1  = (m) -> size(m.aero_model.K8, 1)

    [Derivatives]
    dpos     = (size=3, output=true)
    deuler   = (size=3, output=true)
    ab       = (size=3, output=true)
    dwb      = (size=3, output=true)
    dcirc_w1 = (m) -> size(m.aero_model.K8, 1)

    [Inputs]
    u_cs          = (m) -> length(m.aero_model.controlSurfaces)
    du_cs         = (m) -> length(m.aero_model.controlSurfaces)
    ddu_cs        = (m) -> length(m.aero_model.controlSurfaces)
    thrust        = 3
    thrust_moment = 3
    vd            = 3
    wd            = 3
    ad            = 3
    dwd           = 3
    rs            = (m) -> 3 * m.aero_model.boundMesh.sizes.totalVertices
    vs            = (m) -> 3 * m.aero_model.boundMesh.sizes.totalVertices
    as            = (m) -> 3 * m.aero_model.boundMesh.sizes.totalVertices
    rho           = 1

    [Outputs]
    Fb          = 3
    Mb          = 3
    Fb_unsteady = 3
    Cb          = 6
    Cs          = 6
    Fmp         = (m) -> 6 * length(m.aero_model.monitorPoints)
    Fa          = (m) -> 3 * (TotalSegments(m.aero_model.boundMesh) + TotalPanels(m.aero_model.boundMesh))
end

function RigidBody6DOF(aero_model::M, mass::T, inertia::SMatrix{3, 3, T, 9}) where {M, T}
    inertia_inv = inv(inertia)
    massAero = SMatrix{6, 6, T, 36}(AddedMass(aero_model))
    return RigidBody6DOF{M, T}(aero_model, mass, inertia, inertia_inv, massAero)
end

# Alias for backwards compatibility
const UnsteadyFlightDynamics = RigidBody6DOF

function AeroPanels.AllocateFMUCaches(model::RigidBody6DOF{M, T}, ::Type{T}) where {M, T}
    return AllocateFMUCaches(model.aero_model, T)
end

function AeroPanels.InitializeFMU!(array::AbstractFMUArray, cache, model::RigidBody6DOF, t::Float64; start_from_trim::Bool=false)
    InitializeFMU!(array, cache, model.aero_model, t; start_from_trim=start_from_trim)
    return nothing
end

function AeroPanels.InitializeSteadyState!(array::AbstractFMUArray, cache, model::RigidBody6DOF, t::Float64)
    InitializeSteadyState!(array, cache, model.aero_model, t)
    return nothing
end

function AeroPanels.EvaluateDerivatives!(array::AbstractFMUArray, cache, model::RigidBody6DOF, t::Float64)
    EvaluateDerivatives!(array, cache, model.aero_model, t)
    
    # 1. Zero out acceleration inputs to compute purely circulation forces F_circ
    fill!(array.ab, 0.0)
    fill!(array.dwb, 0.0)
    EvaluateOutputs!(array, cache, model.aero_model, t)
    
    Fb_circ = SVector{3}(array.Fb) + SVector{3}(array.thrust)
    Mb_circ = SVector{3}(array.Mb) + SVector{3}(array.thrust_moment)
    
    # 2. Solve coupled 6DOF equations of motion using (M_rigid + rho * M_added)
    dp, dv, de, dw = EOM6DOF(SVector{3}(array.vb), SVector{3}(array.wb), SVector{3}(array.euler), Fb_circ, Mb_circ, model, array.rho[1])
    
    array.dpos .= dp
    array.deuler .= de
    array.ab .= dv
    array.dwb .= dw
    
    # 3. Analytical superposition of acceleration added-mass forces onto body forces
    accel_6dof = SA[dv[1], dv[2], dv[3], dw[1], dw[2], dw[3]]
    F_accel_6dof = - array.rho[1] * (model.massAero * accel_6dof)
    array.Fb .+= F_accel_6dof[1:3]
    array.Mb .+= F_accel_6dof[4:6]
    
    return nothing
end

function AeroPanels.EvaluateOutputs!(array::AbstractFMUArray, cache, model::RigidBody6DOF, t::Float64)
    # Evaluates complete aerodynamic outputs (including acceleration added mass forces) using array.ab and array.dwb
    EvaluateOutputs!(array, cache, model.aero_model, t)
    return nothing
end

function AeroPanels.ManualDirectionalDerivative(
    model::RigidBody6DOF{M, T},
    array::AbstractFMUArray,
    outputRefs::AbstractVector{<:Integer},
    inputRefs::AbstractVector{<:Integer},
    inputSeed::AbstractVector;
    cache=nothing
) where {M, T}
    return ManualDirectionalDerivative(model.aero_model, array, outputRefs, inputRefs, inputSeed; cache=cache)
end

"""
    EOM6DOF(vel, omega, euler, F_aero, M_aero, model::RigidBody6DOF{M, T}, rho) where {M, T}

Evaluates the 6DOF rigid body kinematics and dynamics derivatives using the coupled added mass formulation.
"""
function EOM6DOF(vel, omega, euler, F_aero, M_aero, model::RigidBody6DOF{M, T}, rho) where {M, T}
    # 1. Flight Kinematics
    phi, theta, psi = euler[1], euler[2], euler[3]
    s_phi, c_phi = sin(phi), cos(phi)
    s_theta, c_theta = sin(theta), cos(theta)
    s_psi, c_psi = sin(psi), cos(psi)
    
    R_ib = @SMatrix [
        c_theta*c_psi   s_phi*s_theta*c_psi - c_phi*s_psi   c_phi*s_theta*c_psi + s_phi*s_psi;
        c_theta*s_psi   s_phi*s_theta*s_psi + c_phi*c_psi   c_phi*s_theta*s_psi - s_phi*c_psi;
        -s_theta        s_phi*c_theta                       c_phi*c_theta
    ]
    
    dpos = R_ib * vel
    
    # Euler rates kinematic matrix
    sec_theta = abs(c_theta) > 1e-6 ? sec(theta) : sign(c_theta)*1e6
    tan_theta = abs(c_theta) > 1e-6 ? tan(theta) : sign(c_theta)*1e6
    
    deuler = SA[
        omega[1] + s_phi*tan_theta*omega[2] + c_phi*tan_theta*omega[3],
        c_phi*omega[2] - s_phi*omega[3],
        s_phi*sec_theta*omega[2] + c_phi*sec_theta*omega[3]
    ]
    
    # 2. Flight Kinetics
    # Gravity in body frame
    F_grav = model.mass * 9.81 * SA[-s_theta, s_phi*c_theta, c_phi*c_theta]
    
    # Construct 6x6 rigid mass matrix
    mass = model.mass
    inertia = model.inertia
    M_rigid = @SMatrix [
        mass  0.0   0.0   0.0           0.0           0.0;
        0.0   mass  0.0   0.0           0.0           0.0;
        0.0   0.0   mass  0.0           0.0           0.0;
        0.0   0.0   0.0   inertia[1,1]  inertia[1,2]  inertia[1,3];
        0.0   0.0   0.0   inertia[2,1]  inertia[2,2]  inertia[2,3];
        0.0   0.0   0.0   inertia[3,1]  inertia[3,2]  inertia[3,3]
    ]
    
    # Total mass matrix (rigid + scaled added mass)
    M_total = M_rigid + rho * model.massAero
    
    # Coriolis/kinematic forces and moments
    F_kin = model.mass * cross(omega, vel)
    
    M_kin = cross(omega, inertia * omega)
    
    b_RHS = SA[
        F_grav[1] + F_aero[1] - F_kin[1],
        F_grav[2] + F_aero[2] - F_kin[2],
        F_grav[3] + F_aero[3] - F_kin[3],
        M_aero[1] - M_kin[1],
        M_aero[2] - M_kin[2],
        M_aero[3] - M_kin[3]
    ]
    
    # Solve coupled equations of motion
    a_6DOF = M_total \ b_RHS
    
    dvel = SA[a_6DOF[1], a_6DOF[2], a_6DOF[3]]
    domega = SA[a_6DOF[4], a_6DOF[5], a_6DOF[6]]
    
    return dpos, dvel, deuler, domega
end
