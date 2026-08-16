# --- Sample Rectangular Wing Models ---

"""
    RectangularWingSteady(; nc=10, ns=20, chord=1.0, span=6.0)

Constructs a steady-state mirrored rectangular wing model of type `SteadyAeroModel2D{Float64}`.
"""
function RectangularWingSteady(; nc=10, ns=20, chord=1.0, span=6.0)
    surf = AeroSurface2D(nc, ns; chord=(chord, chord), span=span/2.0, sweep=0.0, dihedral=0.0)
    surf_mirrored = Mirror(surf, 2)
    props = AeroModelProperties(c=chord, b=span, S=span * chord)
    
    boundMesh = ThinAeroMesh([surf, surf_mirrored])
    wakeMesh = FixedWakeMesh(boundMesh.ringMesh, boundMesh.sizes, props)
    
    return SteadyAeroModel2D(boundMesh, wakeMesh, props)
end

"""
    RectangularWingUnsteady(; nc=10, ns=20, chord=1.0, span=6.0, nWake=80, wakeLength=20.0)

Constructs an unsteady mirrored rectangular wing model of type `UnsteadyAeroModel2D{Float64}`.
"""
function RectangularWingUnsteady(; nc=10, ns=20, chord=1.0, span=6.0, nWake=80, wakeLength=20.0)
    surf = AeroSurface2D(nc, ns; chord=(chord, chord), span=span/2.0, sweep=0.0, dihedral=0.0)
    surf_mirrored = Mirror(surf, 2)
    props = AeroModelProperties(c=chord, b=span, S=span * chord)
    
    boundMesh = ThinAeroMesh([surf, surf_mirrored])
    wakeMesh = FixedWakeMesh(boundMesh.ringMesh, boundMesh.sizes, props; nWake=nWake, wakeLength=wakeLength)
    
    return UnsteadyAeroModel2D(boundMesh, wakeMesh, props, 1.0)
end
