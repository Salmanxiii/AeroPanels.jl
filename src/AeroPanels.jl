module AeroPanels

using LinearAlgebra
using StaticArrays
using GeometryBasics
using SparseArrays
using Polyester
using DocStringExtensions

import Base: show, size, iterate, length, @kwdef

include("Indexing.jl")
include("PanelProperties.jl")
include("Segments.jl")
include("AeroSurface2D.jl")
include("Meshing.jl")
include("Influence.jl")
include("Steady/SteadyAeroModel2D.jl")
include("Steady/solve.jl")
include("ControlSurfaces.jl")
include("MonitorPoints.jl")
include("Unsteady/UnsteadyAeroModel2D.jl")
include("Misc.jl")


# Exports
export Sizes, IndicesMatrix, SelectionOperator
export AeroModelProperties, FlowAxis, PanelProperties, SegmentProperties, ProcessSegments
export SegmentCirculation, SegmentInducedVelocity, SegmentForce, NormalWash
export AeroSurface, AeroSurface2D, Mirror, NoSegments
export WakeModel, CreateAeroMesh, RingMesh, FlatWakeMesh
export VORTXL, VORING, Influence, SteadyWakeInfluence
export AeroModel, AeroModel2D, SteadySolution, AICSolve, AerodynamicForces, AeroSolve
export BodyAccelerations, AerodynamicAnglesDerivatives
export UnsteadyAeroModel2D, UnsteadyWakeInfluence, FullWakeFromTransportWakeOperator
export GetFullWakeVector, SolveSteadyCirculation, SolveCirculation, UnsteadyPanelForces
export SolveForces, NumberOfStates
export PlotModel

export GetTotalForces, GetCoefficients, AeroSolve!

function PlotModel end

end
