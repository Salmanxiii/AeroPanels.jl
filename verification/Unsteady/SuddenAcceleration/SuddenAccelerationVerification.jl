using AeroPanels
using LinearAlgebra
using OrdinaryDiffEq
using CairoMakie
using StaticArrays

# --- Case Definition ---
c = 1.0          # Chord [m]
V = 1.0          # Velocity [m/s]
α = 2.0          # Angle of attack [deg]
ρ = 1.225        # Air density [kg/m^3]
nc = 20          # Number of chordwise panels
ns = 1           # 1 spanwise panel (2D)
ar = 1000        # Large aspect ratio for 2D approximation

# Create Model
surf = AeroSurface2D(nc, ns, chord=(c, c), span=ar*c)
props = AeroModelProperties(c=c, b=ar*c, S=c*ar*c, symmXZ=false)
model = UnsteadyAeroModel2D([surf], props, V, nWake=100, wakeLength=40.0)

# Sudden Acceleration setup
vb = SVector(V*cos(deg2rad(α)), 0.0, V*sin(deg2rad(α)))
dvb = SVector(0.0, 0.0, 0.0) # Constant velocity
b = AeroPanels.NormalWash(vb, model)

# ODE setup
u0 = zeros(NumberOfStates(model))
tspan = (0.0, 50.0 * (c/V)) # 50 chord lengths

function uvlm_dynamics!(du, u, p, t)
    du .= SolveCirculation(u, b, model)
end

prob = ODEProblem(uvlm_dynamics!, u0, tspan)
sol = solve(prob, RK4(), reltol=1e-6, abstol=1e-6)

# Post-processing: Calculate lift coefficient
CL = Float64[]
time_chord = sol.t ./ (c/V) # Non-dimensional time (t*V/c)
s = 2.0 .* time_chord       # Wagner's s (2*V*t/c)

for i in 1:length(sol.t)
    Γw = sol.u[i]
    Fa, Fu = SolveForces(Γw, vb, dvb, model, ρ)
    
    # Total vertical force
    Fz = sum(f[3] for f in Fa) + sum(f[3] for f in Fu)
    
    # Lift coefficient CL = L / (0.5 * ρ * V^2 * S)
    cl = Fz / (0.5 * ρ * V^2 * props.S)
    push!(CL, cl)
end

# Analytical Solution: Wagner's Function (Jones Approximation)
# w(s) = 1 - 0.165*exp(-0.0455*s) - 0.335*exp(-0.3*s)
wagner(s) = 1.0 - 0.165*exp(-0.0455*s) - 0.335*exp(-0.3*s)
CL_steady = 2 * π * deg2rad(α) # Thin airfoil theory steady lift
CL_analytical = CL_steady .* wagner.(s)

# Plotting
fig = Figure(resolution = (800, 600))
ax = Axis(fig[1, 1],
    title = "Sudden Acceleration (Wagner's Problem)",
    xlabel = L"t V / c",
    ylabel = L"C_L")

lines!(ax, time_chord, CL/CL_steady, label="UVLM (AeroPanels.jl)", linewidth=2)
lines!(ax, time_chord, CL_analytical/CL_steady, label="Wagner (Jones approx.)", linestyle=:dash, linewidth=2)

axislegend(ax, position = :rb)

save("Unsteady/SuddenAcceleration/SuddenAcceleration.png", fig)
println("Verification plot saved to verification/Unsteady/SuddenAcceleration/SuddenAcceleration.png")

# Final check
println("Final CL (UVLM): ", CL[end])
println("Final CL (Steady Theory): ", CL_steady)
println("Error: ", abs(CL[end] - CL_steady)/CL_steady * 100, "%")
