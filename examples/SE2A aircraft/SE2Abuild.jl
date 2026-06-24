using AeroPanels
using StaticArrays
using JuliaC
using XML
using CodeTracking

# --- 1. Model Definition ---
# This function defines the aerodynamic model that we want to export as an FMU.
# It MUST return an UnsteadyAeroModel2D.
function create_se2a_model()
    # Wing Geometry
    wingSweep = deg2rad(13.06)
    wingOrigin = [16.425, 0, -1.4]
    wingSpan = [0., 2.04, 3.87, 5.50, 6.96, 8.33, 9.63, 10.85, 12.00,
                13.08, 14.11, 15.07, 15.98, 16.84, 17.65, 18.41, 19.13, 19.81, 20.44, 21.04, 21.61]
    wingChord = [6.3994, 5.7177, 5.1084, 4.5648, 4.1853, 3.9475, 3.7231, 3.5115, 3.3119, 3.1237, 
                 2.9461, 2.7787, 2.6186, 2.4673, 2.3248, 2.1906, 2.0641, 1.9449, 1.8325, 1.7282, 1.6300]
    xDivsWing = [0.0, 0.066, .133, 0.2, 0.3, 0.4, .5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1.0]


    stbdWing = AeroSurface2D(length(xDivsWing)-1, length(wingSpan)-1;
                                chord=wingChord, span=wingSpan,
                                xDivs=xDivsWing, sweep=wingSweep,
                                dihedral=0, origin=wingOrigin)
    fullWing = MakeSymmetricY(stbdWing)

    # HTail Geometry
    htSweep = 0.4566
    htOrigin = [32.5, 0, 1.25]
    htChord = [3.7500, 3.3447, 2.9828, 2.6601, 2.3723, 2.1156, 1.8867, 1.6826, 1.5000]
    htSpan = [0.0000, 1.1694, 2.2124, 3.1425, 3.9719, 4.7117, 5.3714, 5.9597, 6.4843]
    xDivsTail = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.775, 0.83125, 0.8875, 0.94375, 1.0]
    stbdTail = AeroSurface2D(length(xDivsTail)-1, length(htSpan)-1;
                                chord=htChord, span=htSpan,
                                xDivs=xDivsTail, sweep=htSweep,
                                dihedral=0, origin=htOrigin)
    fullTail = MakeSymmetricY(stbdTail)

    # VTail Geometry
    vtOrigin = [30, 0., 1.25]
    vtSweep = [0.484, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537]
    vtChord = [6.00, 5.23, 4.55, 3.95, 3.43, 2.98, 2.59, 2.25, 1.95]
    vtSpan = [0.00, 1.52, 2.83, 3.98, 4.97, 5.83, 6.58, 7.23, 7.80]
    VTail = AeroSurface2D(length(xDivsTail)-1, length(vtSpan)-1;
                                chord=vtChord, span=vtSpan,
                                xDivs=xDivsTail, sweep=vtSweep,
                                dihedral=pi/2, origin=vtOrigin)

    aeroSurfaces = [fullWing, fullTail, VTail]

    # Control Surfaces
    pail   = CSDefinition("PortAileron", 1, (10,14), (3,10))
    pflap2 = CSDefinition("PortFlap2", 1, (10,14), (11,17))
    pflap1 = CSDefinition("PortFlap1", 1, (10,14), (19,20))
    sflap1 = CSDefinition("StbdFlap1", 1, (10,14), (22,23))
    sflap2 = CSDefinition("StbdFlap2", 1, (10,14), (25,31))
    sail   = CSDefinition("StbdAileron", 1, (10,14), (32,39))

    pele = CSDefinition("PortElevator", 2, (9,13), (2, 8))
    sele = CSDefinition("StbdElevator", 2, (9,13), (10,16))
    rudd = CSDefinition("Rudder", 3, (9,13), (2,9))

    CSvec = [pail, pflap2, pflap1, sflap1, sflap2, sail, pele, sele, rudd]

    # Load Monitor Points
    wingMPorigin = wingOrigin + SA[wingChord[1]/4, 0, 0] # at 25% chord
    mp1 = MPDefinition("StbdWingRoot", 1, wingMPorigin, (1,14), (1,20))
    mp2 = MPDefinition("PortWingRoot", 1, wingMPorigin, (1,14), (21,41))
    mpVec = [mp1, mp2]

    # Global Properties of AeroModel
    props = AeroModelProperties(
        c=3.668, b=43.2228, S=158.5356, CG=SA[20.13, 0., -0.1018]);

    model = UnsteadyAeroModel2D(aeroSurfaces, props, 1., nWake=80, wakeLength=20.,
     controlSurfaces=CSvec, monitorPoints=mpVec)
    return model
end

fmu_name = "SE2AAircraftFMU"
fmu_dir = @__DIR__
fmu_path = joinpath(fmu_dir, fmu_name * ".fmu")

# --- 2. Build the FMU ---
println("Building FMU for SE2A Aircraft model...")
BuildFMU(create_se2a_model, fmu_dir, fmu_name=fmu_name, clean=true)

println("--- SUCCESS ---")
println("FMU created at: $fmu_path")
