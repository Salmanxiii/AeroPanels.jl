using AeroPanels
using GLMakie
using StaticArrays
using OrdinaryDiffEqTsit5

# Ensure assets directory exists
mkpath(joinpath(@__DIR__, "src", "assets"))

# --- 1. Steady Swept Wing Images ---
nc, ns = 10, 20
chord = (1.0, 1.0)
span = 10.0
sweep_angle = deg2rad(30)
surf = AeroSurface2D(nc, ns, chord=chord, span=span, sweep=sweep_angle)
surf2 = Mirror(surf, 2)
props = AeroModelProperties(c=1.0, b=span, S=span*1.0)
boundMesh = ThinAeroMesh([surf, surf2])
wakeMesh = FixedWakeMesh(boundMesh.ringMesh, boundMesh.sizes, props)
model = SteadyAeroModel2D(boundMesh, wakeMesh, props)

V, α = 10.0, 5.0
vb = AeroPanels.BodyVelocity(V, deg2rad(α))
cache = AeroSolve(vb, model)

# Circulation
fig3 = PlotModel(model, plotWake=true, Γp=cache.Γp, Γw=cache.Γw)
save(joinpath(@__DIR__, "src", "assets", "circulation.png"), fig3)
