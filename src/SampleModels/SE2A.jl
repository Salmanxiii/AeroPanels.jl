# --- Sample SE2A Aircraft Models ---

#   SE2AGeometry()
#
# Internal helper function that constructs and returns the surfaces and global properties
function SE2AGeometry()
    # Wing Geometry
    wingSweep = deg2rad(13.06)
    wingOrigin = [16.425, 0.0, -1.4]
    wingSpan = [0.0, 2.04, 3.87, 5.50, 6.96, 8.33, 9.63, 10.85, 12.00,
                13.08, 14.11, 15.07, 15.98, 16.84, 17.65, 18.41, 19.13, 19.81, 20.44, 21.04, 21.61]
    wingChord = [6.3994, 5.7177, 5.1084, 4.5648, 4.1853, 3.9475, 3.7231, 3.5115, 3.3119, 3.1237, 
                 2.9461, 2.7787, 2.6186, 2.4673, 2.3248, 2.1906, 2.0641, 1.9449, 1.8325, 1.7282, 1.6300]
    xDivsWing = [0.0, 0.066, 0.133, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1.0]

    stbdWing = AeroSurface2D(length(xDivsWing)-1, length(wingSpan)-1;
                                chord=wingChord, span=wingSpan,
                                xDivs=xDivsWing, sweep=wingSweep,
                                dihedral=0.0, origin=wingOrigin)
    fullWing = MakeSymmetricY(stbdWing)

    # HTail Geometry
    htSweep = 0.4566
    htOrigin = [32.5, 0.0, 1.25]
    htChord = [3.7500, 3.3447, 2.9828, 2.6601, 2.3723, 2.1156, 1.8867, 1.6826, 1.5000]
    htSpan = [0.0000, 1.1694, 2.2124, 3.1425, 3.9719, 4.7117, 5.3714, 5.9597, 6.4843]
    xDivsTail = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.775, 0.83125, 0.8875, 0.94375, 1.0]
    stbdTail = AeroSurface2D(length(xDivsTail)-1, length(htSpan)-1;
                                chord=htChord, span=htSpan,
                                xDivs=xDivsTail, sweep=htSweep,
                                dihedral=0.0, origin=htOrigin)
    fullTail = MakeSymmetricY(stbdTail)

    # VTail Geometry
    vtOrigin = [30.0, 0.0, 1.25]
    vtSweep = [0.484, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537]
    vtChord = [6.00, 5.23, 4.55, 3.95, 3.43, 2.98, 2.59, 2.25, 1.95]
    vtSpan = [0.00, 1.52, 2.83, 3.98, 4.97, 5.83, 6.58, 7.23, 7.80]
    VTail = AeroSurface2D(length(xDivsTail)-1, length(vtSpan)-1;
                                chord=vtChord, span=vtSpan,
                                xDivs=xDivsTail, sweep=vtSweep,
                                dihedral=pi/2, origin=vtOrigin)

    aeroSurfaces = [fullWing, fullTail, VTail]

    # Global Properties of AeroModel
    props = AeroModelProperties(
        c=3.668, b=43.2228, S=158.5356, rRef=SA[20.13, 0.0, -0.1018])

    return (aeroSurfaces, props, wingOrigin, wingChord[1])
end

function SE2AFeatures(boundMesh, wingOrigin, rootChord)
    # Control Surfaces
    pail   = CreateControlSurface(boundMesh, "PortAileron", surfaceIdx=1, nc=(10,13), ns=(3,10))
    pflap2 = CreateControlSurface(boundMesh, "PortFlap2", surfaceIdx=1, nc=(10,13), ns=(11,17))
    pflap1 = CreateControlSurface(boundMesh, "PortFlap1", surfaceIdx=1, nc=(10,13), ns=(19,20))
    sflap1 = CreateControlSurface(boundMesh, "StbdFlap1", surfaceIdx=1, nc=(10,13), ns=(22,23))
    sflap2 = CreateControlSurface(boundMesh, "StbdFlap2", surfaceIdx=1, nc=(10,13), ns=(25,31))
    sail   = CreateControlSurface(boundMesh, "StbdAileron", surfaceIdx=1, nc=(10,13), ns=(32,39))

    pele = CreateControlSurface(boundMesh, "PortElevator", surfaceIdx=2, nc=(9,12), ns=(2, 8))
    sele = CreateControlSurface(boundMesh, "StbdElevator", surfaceIdx=2, nc=(9,12), ns=(10,16))
    rudd = CreateControlSurface(boundMesh, "Rudder", surfaceIdx=3, nc=(9,12), ns=(2,8))

    CSvec = [pail, pflap2, pflap1, sflap1, sflap2, sail, pele, sele, rudd]

    # Load Monitor Points
    wingMPorigin = wingOrigin + SA[rootChord/4, 0.0, 0.0] # at 25% chord
    mp1 = CreateMonitorPoint(boundMesh, "StbdWingRoot", surfaceIdx=1, nc=(1,14), ns=(1,20), origin=wingMPorigin)
    mp2 = CreateMonitorPoint(boundMesh, "PortWingRoot", surfaceIdx=1, nc=(1,14), ns=(21,41), origin=wingMPorigin)
    mpVec = [mp1, mp2]

    return CSvec, mpVec
end

"""
    SE2ASteady()

Constructs the canonical steady-state SE2A aircraft model of type `SteadyAeroModel2D{Float64}`.
"""
function SE2ASteady()
    aeroSurfaces, props, wingOrigin, rootChord = SE2AGeometry()
    boundMesh = ThinAeroMesh(aeroSurfaces)
    wakeMesh = FixedWakeMesh(boundMesh.ringMesh, boundMesh.sizes, props)
    CSvec, mpVec = SE2AFeatures(boundMesh, wingOrigin, rootChord)
    
    return SteadyAeroModel2D(boundMesh, wakeMesh, props; controlSurfaces=CSvec, monitorPoints=mpVec)
end

"""
    SE2AUnsteady(; nWake=80, wakeLength=20.0)

Constructs the canonical unsteady SE2A aircraft model of type `UnsteadyAeroModel2D{Float64}`.
"""
function SE2AUnsteady(; nWake=80, wakeLength=20.0, wakeExpansionR=1.0)
    aeroSurfaces, props, wingOrigin, rootChord = SE2AGeometry()
    boundMesh = ThinAeroMesh(aeroSurfaces)
    wakeMesh = FixedWakeMesh(boundMesh.ringMesh, boundMesh.sizes, props; nWake=nWake, wakeLength=wakeLength, wakeExpansionR=wakeExpansionR)
    CSvec, mpVec = SE2AFeatures(boundMesh, wingOrigin, rootChord)
    
    return UnsteadyAeroModel2D(boundMesh, wakeMesh, props, 1.0; controlSurfaces=CSvec, monitorPoints=mpVec)
end

"""
    SE2AUnsteadyFlight(; nWake=80, wakeLength=20.0)

Constructs the coupled unsteady flight dynamics model (`UnsteadyFlightDynamics`) for the SE2A aircraft.
"""
function SE2AUnsteadyFlight(; nWake=80, wakeLength=20.0, wakeExpansionR=1.0)
    aero_model = SE2AUnsteady(; nWake=nWake, wakeLength=wakeLength, wakeExpansionR=wakeExpansionR)
    mass = 64158.0
    inertia = @SMatrix [
         1.191854814680426e+06    2.745065883456391    -9.435430560624978e+04;
         2.745065883456391        3.392996536125001e+06    -0.843080668520315;
        -9.435430560624978e+04    -0.843080668520315         4.499840235119487e+06
    ]
    return UnsteadyFlightDynamics(aero_model, mass, inertia)
end
