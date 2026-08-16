module AeroPanels

using LinearAlgebra
using StaticArrays
using GeometryBasics
using SparseArrays
using Polyester
using DocStringExtensions
using OrdinaryDiffEqTsit5: ODEProblem, solve
using DiffEqCallbacks
using ForwardDiff

import Base: show, size, iterate, length, @kwdef

abstract type AeroModel end

include("Interfaces/FMUInterfaces.jl")
include("Interfaces/Simulation.jl")

include("Indexing.jl")
include("Geometry/AbstractAeroMesh.jl")
include("Geometry/AeroSurface2D.jl")
include("Geometry/ThinAeroMesh.jl")
include("Geometry/Segments.jl")
include("Geometry/Meshing.jl")

GetBoundMesh(model::AeroModel) = GetBoundMesh(model.boundMesh)
GetWakeMesh(model::AeroModel) = GetWakeMesh(model.wakeMesh)
GetWakeMesh(mesh) = mesh
include("Influence.jl")
include("Geometry/ControlSurfaces.jl")
include("Geometry/MonitorPoints.jl")
include("Misc.jl")
include("Models/SteadyVLM.jl")
include("Models/SteadyVLMSolve.jl")
include("Models/UnsteadyVLM.jl")
include("Models/UnsteadyVLMSolve.jl")
include("Models/Dynamics/RigidBody6DOF.jl")

include("SampleModels/RectangularWing.jl")
include("SampleModels/SE2A.jl")

function BuildFMU end

# Exports
export Sizes, IndicesMatrix, SelectionOperator
export AeroModelProperties, FlowAxis, SegmentProperties, ProcessSegments
export AbstractAeroMesh, ThinAeroMesh, FixedWakeMesh
export GetCollocationPoints, GetMidpoints, GetNormalVectors, GetPanelAreas
export GetBoundPanelCount, GetBoundMesh, GetWakeMesh, GetControlSurfaces, GetMonitorPoints
export SegmentCirculation, SegmentInducedVelocity, SegmentForce, NormalWash
export AeroSurface, AeroSurface2D, Mirror, MakeSymmetricY, TotalSegments, TotalPanels
export CreateAeroMesh, RingMesh, FixedWakeMesh
export VORTXL, VORING, Influence
export AeroModel, SteadyAeroModel2D, AeroSolve, RigidBody6DOF, UnsteadyFlightDynamics, AddedMass, EOM6DOF
export BodyAccelerations, AerodynamicAnglesDerivatives, BodyVelocity
export UnsteadyAeroModel2D, UnsteadyWakeInfluence, FullWakeFromTransportWakeOperator
export GetFullWakeVector, SolveSteadyCirculation, SolveCirculation, UnsteadyPanelForces
export SolveForces, NumberOfStates, CreateControlSurface, ControlSurface, CreateMonitorPoint, MonitorPoint
export PlotModel, GetStabilityCoefficients
export RectangularWingSteady, RectangularWingUnsteady, SE2ASteady, SE2AUnsteady, SE2AUnsteadyFlight

export GetTotalForces, GetStabilityDerivatives, AeroSolve!, CreateCacheArrays, MonitorPointLoads!
export BuildFMU, UnsteadyAeroCache, InputCache, OutputCache, StateCache
export AbstractFMUModel, AbstractFMUArray, CreateFMUArray, GetTotalSize, @fmu_model, FMULayout, GetFMULayout, EvaluateDerivatives!, EvaluateOutputs!, InitializeFMU!, InitializeSteadyState!, AllocateFMUCaches, simulate, GetDirectionalDerivative!, GetFMISparsity

function PlotModel end

end
