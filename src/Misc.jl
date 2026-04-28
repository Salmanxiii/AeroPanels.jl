GeometryToBodyAxis(Vec, props::AeroModelProperties) = props.bodyFixedCS * Vec

BodyFixedToStabilityAxis(Vec, α) = Point3(Vec[1] * cos(α) + Vec[3] * sin(α),
 Vec[2], -Vec[1] * sin(α) + Vec[3] * cos(α))

BodyVelocity(V, α, β=0.) = Point3(V * cos(α) * cos(β), V * sin(β), V * sin(α) * cos(β))

function BodyAccelerations(Vt, dVt, α, dα, β=0., dβ=0.)
    sα, cα = sincos(α)
    sβ, cβ = sincos(β)

    dU = dVt*cα*cβ - Vt*dα*sα*cβ - Vt*dβ*cα*sβ
    dv = dVt*sβ    + Vt*dβ*cβ
    dW = dVt*sα*cβ + Vt*dα*cα*cβ - Vt*dβ*sα*sβ
    return Point3(dU, dv, dW)
end

AerodynamicAngles(vb) = (atan(vb[3], vb[1]), asin(vb[2]/norm(vb)), norm(vb))
AerodynamicAnglesDegree(vb) = (atand(vb[3], vb[1]), asind(vb[2]/norm(vb)), norm(vb))

function AerodynamicAnglesDerivatives(vb, dvb)
    U, v, W = vb
    U̇, v̇, Ẇ = dvb
    U2W2 = U^2 + W^2
    Vt = norm(vb)

    dVt = dot(V, V̇) / Vt
    dα = (U * Ẇ - W * U̇) / U2W2
    dβ = (v̇ * Vt - v * dVt) / (Vt * sqrt(U2W2))
    return dα, dβ, dVt
end

GeometryToStabilityAxis(Vec, α, props::AeroModelProperties) =  BodyFixedToStabilityAxis(GeometryToBodyAxis(Vec, props), α)
