abstract type AeroSurface{T} end

"""
$(TYPEDEF)

A 2D aerodynamic surface defined by a mesh of points in 3D space.

$(TYPEDFIELDS)
"""
struct AeroSurface2D{T<:Real} <: AeroSurface{T}
    X::Matrix{T}
    Y::Matrix{T}
    Z::Matrix{T}
    nc::Int
    ns::Int
end

"""
    AeroSurface2D(nc, ns; chord, span, xDivs, camber, sweep, dihedral, origin)

High-level constructor for an `AeroSurface2D`.
"""
function AeroSurface2D(nc::Int, ns::Int;
   chord::Union{Real, Tuple{Real,Real}, AbstractVector} = 1,
   span::Union{Real, AbstractVector} = 1,
   xDivs::AbstractVector = LinRange(0, 1, nc+1),
   camber::Union{AbstractArray, AbstractMatrix} = zeros(nc+1),
   sweep::Union{Real, AbstractVector} = 0,
   dihedral::Union{Real, AbstractVector} = 0,
   origin::AbstractVector = SA[0,0,0])

   if typeof(span) <: Real
       span = LinRange(0, span, ns+1)
   end

   if typeof(chord) <: Real
       chord = repeat([chord], ns+1)
   elseif typeof(chord) <: Tuple{Real,Real}
       cRoot = chord[1]
       cTip = chord[2]
       m = (cTip - cRoot)/span[end]
       chord = [m * s + cRoot for s in span]
   end

   if typeof(camber) <: AbstractVector
       camber = repeat(collect(camber), 1, ns+1)
   end

   if typeof(sweep) <: Real
       sweep = repeat([sweep], ns)
   end

   if typeof(dihedral) <: Real
       dihedral = repeat([dihedral], ns)
   end

   return AeroSurface2D(nc, ns, chord, span, xDivs, camber, sweep , dihedral, origin)
end

function AeroSurface2D(nc::Int, ns::Int,
    chord::AbstractVector{B}, span::AbstractVector{C},
    xDivs::AbstractVector{A}, camber::AbstractMatrix{D},
    sweep::AbstractVector{E}, dihedral::AbstractVector{F},
    origin::AbstractVector{G}) where {A<:Real, B<:Real,
        C<:Real, D<:Real, E<:Real, F<:Real, G<:Real}

    T = promote_type(A,B,C,D,E,F,G)
    
    (xDivs, chord, span, sweep , dihedral, camber, origin) = 
    (T.(xDivs), T.(chord), T.(span),T.(sweep), T.(dihedral), T.(camber), T.(origin))

    return AeroSurface2D(nc, ns, chord, span, xDivs, camber, sweep , dihedral, origin)
end

function AeroSurface2D(nc::Int, ns::Int,
    chord::AbstractVector{T}, span::AbstractVector{T},
    xDivs::AbstractVector{T}, camber::AbstractMatrix{T},
    sweep::AbstractVector{T}, dihedral::AbstractVector{T},
    origin::AbstractVector{T}) where {T<:Real}

    if length(xDivs) != nc+1; error("argument 'xDivs' should be vector of size nc+1") end
    if length(span) != ns+1; error("argument 'span' should be vector of size ns+1") end
    if length(chord) != ns+1; error("argument 'chord' should be vector of size ns+1") end
    if size(camber) != (nc+1, ns+1); error("argument 'camber' should be matrix of size (nc+1, ns+1)") end
    if length(sweep) != ns; error("argument 'sweep' should be vector of size ns") end
    if length(dihedral) != ns; error("argument 'dihedral' should be vector of size ns") end
    if length(origin) != 3; error("argument 'origin' should be vector of size 3") end

    mSize = (nc+1, ns+1)
    X = zeros(T, mSize) .+ origin[1]
    Y = zeros(T, mSize) .+ origin[2] 
    Z = zeros(T, mSize) .+ origin[3]

    for (j, x) in enumerate(xDivs)
        X[j, 1] = X[j, 1] + x*chord[1]
        Z[j, 1] = Z[j, 1] + camber[j, 1] * chord[1]
    end

    xLE = 0
    for (itemp, s) in enumerate(span[2:end])
        i = itemp + 1
        c = chord[i]
        Λ = sweep[i-1]
        Γ = dihedral[i-1]
        Δs = s - span[i-1]
        xLE = xLE + tan(Λ) * Δs
        for (j, x) in enumerate(xDivs)
            zcamber = camber[j, i] * c
            X[j, i] = X[j, i] + x*c + xLE
            Y[j, i] = Y[j, i-1] + Δs*cos(Γ) - zcamber*sin(Γ)
            #Y[j, i] = Y[j, i-1] + Δs
            Z[j, i] = Z[j, i-1] + Δs*sin(Γ) + zcamber*cos(Γ)
        end
    end
    return AeroSurface2D(X,Y,Z, nc, ns)
end

MirrorX(s::AeroSurface2D) = AeroSurface2D(-s.X, s.Y, s.Z, s.nc, s.ns)
MirrorY(s::AeroSurface2D) = AeroSurface2D(s.X, -s.Y, s.Z, s.nc, s.ns)
MirrorZ(s::AeroSurface2D) = AeroSurface2D(s.X, s.Y, -s.Z, s.nc, s.ns)

"""
    Mirror(surface, axis)

Mirror an `AeroSurface2D` along a given axis (1=X, 2=Y, 3=Z).
"""
function Mirror(surface::AeroSurface2D{T}, axis) where T
    if axis==1
        return MirrorX(surface)
    elseif axis==2
        return MirrorY(surface)
    elseif axis==3
        return MirrorZ(surface)
    else error("invalid value of argument 'axis' ") end
end

function Base.size(surface::AeroSurface2D)
    return (surface.nc, surface.ns)
end

function NoSegments(surface::AeroSurface2D)
    nc, ns = Base.size(surface)
    nSegX = nc * (ns + 1)
    nSegY = ns * (nc + 1)
    return (nSegX, nSegY)
end

function Base.show(io::IO, surface::AeroSurface2D{T}) where {T}
    nc, ns = size(surface)
    print(io,"(nc=",nc, ", ns=", ns, ')')
end
function Base.show(io::IO, ::MIME"text/plain", surf::AeroSurface2D{T}) where {T}
    print(io, "AeroSurface2D{$T}:\n   ", surf)
end
