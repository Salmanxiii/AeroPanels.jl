```@meta
CurrentModule = AeroPanels
```

# AeroPanels.jl

*A fast, continuous-time Unsteady Vortex Lattice Method (UVLM) solver written in Julia.*

Documentation for [AeroPanels.jl](https://github.com/Julius/AeroPanels.jl).

AeroPanels.jl is a potential flow aerodynamic solver designed for high-performance and seamless integration with continuous-time simulations like aeroelastic flutter and flight dynamics. It provides both a traditional steady-state solver and an advanced unsteady state-space formulation.

![Circulation](assets/circulation.png)

## Quick Start: Steady Wing Analysis

Here is a simple example showing how to create a wing, generate its mesh, and solve for the steady aerodynamic forces (Lift and Drag).

```julia
using AeroPanels
using StaticArrays

# 1. Define the Wing Geometry (e.g., a simple swept wing)
nc, ns = 10, 20
chord = (1.0, 1.0)
span = 10.0
sweep_angle = deg2rad(30)
surf = AeroSurface2D(nc, ns, chord=chord, span=span, sweep=sweep_angle)

# 2. Set Model Properties
# Define reference chord (c), reference span (b), and reference area (S)
props = AeroModelProperties(c=1.0, b=span, S=span*1.0)

# 3. Assemble the Aerodynamic Model
model = AeroModel2D([surf], props)

# 4. Define Flight Condition & Solve
V = 10.0   # Freestream speed [m/s]
α = 5.0    # Angle of attack [deg]
ρ = 1.225  # Air density [kg/m^3]

sol = AeroSolve(V, α, model, ρ)

# 5. Extract Results
println("Lift Coefficient (CL): ", round(sol.CL, digits=4))
println("Drag Coefficient (CD): ", round(sol.CD, digits=4))
```

### Next Steps

- Read the [Steady Theory](steady_theory.md) and [Unsteady Theory](unsteady_theory.md) pages to understand the mathematics behind the solvers.
- Check out the [Examples](examples.md) page for more advanced usage, including unsteady Wagner problem simulations and Makie visualizations.
- Browse the [API Reference](api.md) for detailed function documentation.

```@index
```
