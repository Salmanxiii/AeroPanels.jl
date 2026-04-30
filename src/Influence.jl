
################################## Aerodynamics #########################################

function VORTXL(xyzp, xyz1, xyz2, Γ)
    r0 = xyz2 .- xyz1
    r1 = xyzp .- xyz1
    r2 = xyzp .- xyz2
    δ = 1e-10
    r1Norm = norm(r1)
    r2Norm = norm(r2)
    temp1 = sum(cross(r1, r2) .^ 2)
    if r1Norm < δ || r2Norm < δ || temp1 < δ
        q =  SA[0.,0.,0.]
    else
        temp2 =  cross(r1,r2) / temp1
        q = Γ / 4 / π * dot(r0, r1/r1Norm - r2/r2Norm) * temp2
    end
    return q
end

function MirrorY(r)
    return typeof(r)(r[1], -r[2], r[3])
end

function VORTXL(xyzp, xyz1, xyz2, Γ, symmXZ)
    uvw = VORTXL(xyzp, xyz1, xyz2, Γ)
    if symmXZ == true
        uvwMirror = MirrorY(VORTXL(MirrorY(xyzp), xyz1, xyz2, Γ))
        uvw = uvw + uvwMirror
    end
    return uvw
end

function VORING(r0, r1, r2, r3, r4, Γ, symmXZ)
    uvw1 = VORTXL(r0, r1, r2, Γ, symmXZ)
    uvw2 = VORTXL(r0, r2, r3, Γ, symmXZ)
    uvw3 = VORTXL(r0, r3, r4, Γ, symmXZ)
    uvw4 = VORTXL(r0, r4, r1, Γ, symmXZ)
    return uvw1 + uvw2 + uvw3 + uvw4
end

function Influence3!(AIC::Matrix{Point3{T}}, rCollocation::Vector{Point3{T}}, mesh, symmXZ) where T
    n = length(rCollocation)
    m = length(mesh)
    @batch for i in 1:n
        r0 = rCollocation[i]
        for j in 1:m
            r1, r2, r3, r4 = mesh[j]
            AIC[i,j] = VORING(r0, r1, r2, r3, r4, 1.0, symmXZ)
        end
    end
end

function Influence3(rCollocation::Vector{Point3{T}}, mesh, symmXZ) where T
    AIC = Matrix{Point3{T}}(undef, length(rCollocation), m = length(mesh))
    Influence!(AIC, rCollocation, mesh, symmXZ)
    return AIC
end

function Influence(rCollocation::Vector{Point3{T}}, normals::Vector{Point3{T}}, mesh, symmXZ) where T
    n = length(rCollocation)
    m = length(mesh)
    AIC = Matrix{T}(undef, n, m)    
    @batch for i in 1:n
        r0, nrml = rCollocation[i], normals[i]
        for j in 1:m
            r1, r2, r3, r4 = mesh[j]
            v = VORING(r0, r1, r2, r3, r4, 1.0, symmXZ)
            AIC[i,j] = dot(v, nrml)
        end
    end
    return AIC
end