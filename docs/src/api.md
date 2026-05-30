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
