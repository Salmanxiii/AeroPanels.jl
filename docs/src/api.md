# API Reference

```@meta
CurrentModule = AeroPanels
```

## Models

```@docs
AeroModel2D
UnsteadyAeroModel2D
AeroModelProperties
PanelProperties
```

## Solvers

```@docs
AeroSolve
AeroSolve!
PrepareCacheForStates!
PrepareCacheForOutputs!
Outputs!
CreateCacheArrays
GetStabilityCoefficients
GetTotalForces
AddSteadyKinematics!
AddUnsteadyKinematics!
CalculateNormalwash!
CalculateDNormalwash!
SolveCirculation!
SolveUnsteadyCirculation!
CalculateAerodynamicForce!
CalculateUnsteadyAerodynamicForce!
MonitorPointLoads!
```

## Secondary Functions

```@docs
UnsteadyAeroCache
InputCache
OutputCache
StateCache
AeroSurface2D
Mirror
SegmentProperties
SegmentCirculation
SegmentCirculation!
ProcessSegments
NumberOfStates
UnsteadyWakeInfluence
```

## Indexing and Utility

```@docs
Sizes
TEPanelIndex
LEPanelIndex
PanelIndex
```
