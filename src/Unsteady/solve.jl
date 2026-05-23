
# vb: body velocity, ab: acceleration, ωb: angular velocity, dωb: angular acceleration,
# δc: control surface deflection, dδc: CS deflection derivative, ddδc: CS deflection double derivaitve,
# vd: disturbance velocity at mesh vertices, ad: disturbance acceleration
function NormalWash(vb::AbstractVector{A}, ab::AbstractVector{B}, ωb::AbstractVector{C}, dωb::AbstractVector{D},
     δc::AbstractVector{E}, dδc::AbstractVector{F}, ddδc::AbstractVector{G},
     vd::AbstractVector, ad::AbstractVector, model::UnsteadyAeroModel2D)
    
    n = model.sizes.totalPanels
    T = promote_type(A, B, C, D, E, F, G)
    b = zeros(T, n)
    db = copy(b)
    vertices = convert(Vector{Point3{T}}, copy(coordinates(model.mesh)))
    velocities = zeros(Point3{T}, length(vertices))
    accelerations = zeros(Point3{T}, length(vertices))

    NormalWash!(b, db, vertices, velocities, accelerations,
        vb, ab, ωb, dωb, δc, dδc, ddδc, vd, ad, model)

end

function NormalWash!(b, db, vertices, velocities, accelerations,
        vb, ab, ωb, dωb, δc, dδc, ddδc, vd, ad, model)
    
    CSDisplacement!(vertices, δc, model)
    CSVelocity!(velocities, dδc, model)
    CSAcceleration!(accelerations, ddδc, model)

    bodyAxis = model.modelProperties.bodyFixedCS # body axis coordinate
    CG = model.modelProperties.CG
    # velocities in geometry frame
    (vb_g, ab_g, ωb_g, dωb_g) = (bodyAxis' * component for component in (vb, ab, ωb, dωb))

    for i in eachindex(vertices)
        r = vertices[i] - CG
        velocities[i] = velocities[i] + vb_g + vd[i] + cross(ωb_g, r)
        accelerations[i] = accelerations[i] + ab_g + ad[i] + cross(dωb_g, r) # centripetal acceleration ignored
    end

    # Normalwash calculation: calculate new normals and interpolate velocity at collocation point
    for (s, nc, ns) in model.sizes
        @batch for i in 1:nc
            for j in 1:ns
                p = PanelIndex(s, i, j, model.sizes)
                i1, i2, i3, i4 = PanelVertexIndices(s, i, j, model.sizes)

                r1, r2, r3, r4 = vertices[i1], vertices[i2], vertices[i3], vertices[i4]
                v1, v2, v3, v4 = velocities[i1], velocities[i2], velocities[i3], velocities[i4]
                a1, a2, a3, a4 = accelerations[i1], accelerations[i2], accelerations[i3], accelerations[i4]

                normalVec = cross(r3 - r1, r2 - r4)
                normalNorm = norm(normalVec)
                normalUnitVec = normalVec / normalNorm
                v = 0.125*(v1 + v2) + 0.375*(v3 + v4)
                a = 0.125*(a1 + a2) + 0.375*(a3 + a4)
                b[p] = - dot(v, normalUnitVec)

                dNormalVec = cross(v3 - v1, r2 - r4) + cross(r3 - r1, v2 - v4)
                dNormalUnitVec = (dNormalVec - normalUnitVec * dot(normalUnitVec, dNormalVec)) / normalNorm
                db[p] = -(dot(a, normalUnitVec) + dot(v, dNormalUnitVec))
            end
        end
    end
    return b, db
end


function SolveCirculation!(dΓw0, Γb, dΓb, Γw, Γs, Γw0, b, model::UnsteadyAeroModel2D)
    dΓw0 .= model.K8*Γw0 .+ model.K9*b
    Γb .= model.L3*Γw0 .- model.L4*b
    dΓb .= model.L5*Γw0 .+ model.L6*b .- model.L4*db
    Γw .= model.L9 * Γw0 + model.L10 * b
    SegmentCirculation!(Γs, Γb, model.segmentProps)
end

"""
$(SIGNATURES)

Calculate the unsteady contribution to aerodynamic forces based on circulation time derivatives.
"""
function SolveForces(ρ, dΓb, model)
    Fu .= ρ .* dΓb .* model.panelProperties.area .* model.panelProperties.normal
    return Fu
end

function UnsteadyPanelForces(Γw, b, ρ, model)
    dΓb = model.L5*Γw .+ model.L6*b
    Fu = ρ .* dΓb .* model.panelProperties.area .* model.panelProperties.normal
    return Fu, dΓb
end

"""
$(SIGNATURES)

Return total aerodynamic forces (Steady + Unsteady) for given states and kinematics.
"""
function SolveForces(Γw, vb, dvb, model::UnsteadyAeroModel2D, ρ = 1.225)
    b = NormalWash(vb, model)
    db = NormalWash(dvb, model)
    Fu, _ = UnsteadyPanelForces(Γw, b, ρ, model, db)
    Fa, _ = SteadyForce(Γw, b, vb, ρ, model)
    return Fa, Fu
end