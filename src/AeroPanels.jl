module AeroPanels

using LinearAlgebra
using StaticArrays
using GeometryBasics
using SparseArrays
using Polyester

import Base: show, size, iterate, length, @kwdef

include("Indexing.jl")
include("Properties.jl")
include("AeroSurface2D.jl")
include("Meshing.jl")
include("Influence.jl")
include("SteadyAeroModel2D.jl")
include("Misc.jl")

# Exports
export Sizes, IndicesMatrix, SelectionOperator
export AeroModelProperties, FlowAxis, PanelProperties, SegmentProperties, ProcessSegments
export SegmentCirculation, SegmentInducedVelocity, SegmentForce
export AeroSurface, AeroSurface2D, Mirror, NoSegments
export WakeModel, CreateAeroMesh, RingMesh, FlatWakeMesh
export VORTXL, VORING, Influence, SteadyWakeInfluence
export AeroModel, AeroModel2D, SteadySolution, AICSolve, AerodynamicForces, AeroSolve
export BodyAccelerations, AerodynamicAnglesDerivatives

end
