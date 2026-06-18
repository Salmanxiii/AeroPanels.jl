module AeroPanels

using LinearAlgebra
using StaticArrays
using GeometryBasics
using SparseArrays
using Polyester
using DocStringExtensions
using OrdinaryDiffEqTsit5: ODEProblem, solve
using ForwardDiff

import Base: show, size, iterate, length, @kwdef

abstract type AeroModel end

include("Indexing.jl")
include("PanelProperties.jl")
include("Segments.jl")
include("AeroSurface2D.jl")
include("Meshing.jl")
include("Influence.jl")
include("ControlSurfaces.jl")
include("MonitorPoints.jl")
include("Steady/SteadyAeroModel2D.jl")
include("Steady/solve.jl")
include("Unsteady/UnsteadyAeroModel2D.jl")
include("Unsteady/solve.jl")
include("Misc.jl")

function BuildFMU end

# Exports
export Sizes, IndicesMatrix, SelectionOperator
export AeroModelProperties, FlowAxis, PanelProperties, SegmentProperties, ProcessSegments
export SegmentCirculation, SegmentInducedVelocity, SegmentForce, NormalWash
export AeroSurface, AeroSurface2D, Mirror, MakeSymmetricY, NoSegments
export WakeModel, CreateAeroMesh, RingMesh, FlatWakeMesh
export VORTXL, VORING, Influence, SteadyWakeInfluence
export AeroModel, AeroModel2D, AeroSolve
export BodyAccelerations, AerodynamicAnglesDerivatives, BodyVelocity
export UnsteadyAeroModel2D, UnsteadyWakeInfluence, FullWakeFromTransportWakeOperator
export GetFullWakeVector, SolveSteadyCirculation, SolveCirculation, UnsteadyPanelForces
export SolveForces, NumberOfStates, CSDefinition, ControlSurface, MPDefinition
export PlotModel, GetStabilityCoefficients

export GetTotalForces, GetStabilityDerivatives, AeroSolve!, CreateCacheArrays, MonitorPointLoads!
export AeroInputs, BuildFMU, UnsteadyAeroCache

function PlotModel end

end
