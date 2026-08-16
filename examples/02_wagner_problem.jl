# # Unsteady Simulation (Wagner Problem)
#
# This example simulates the lift buildup on a flat plate following a sudden acceleration 
# (impulsive start), comparing the numerical result against the classical Wagner problem.

using AeroPanels
using OrdinaryDiffEqTsit5
using StaticArrays
using GLMakie

# # 1. Simulation Setup
# We define a high-aspect ratio wing to approximate 2D behavior and setup the dynamics.
c, V, α = 1.0, 1.0, 2.0
surf = AeroSurface2D(20, 1, chord=(c, c), span=100.0)
props = AeroModelProperties(c=c, b=100.0, S=100.0)
boundMesh = ThinAeroMesh([surf])
wakeMesh = FixedWakeMesh(boundMesh.ringMesh, boundMesh.sizes, props; nWake=500, wakeLength=20.0)
model = UnsteadyAeroModel2D(boundMesh, wakeMesh, props, V)
tspan = (0.0, 20.0)

# # 2. Input Function
# The input function modifies the AeroInputs struct in-place for each time step.
# For the Wagner problem, the velocity is constant after the impulsive start.
function wagner_inputs!(inputs, t)
    inputs.vb .= SVector(V*cos(deg2rad(α)), 0.0, V*sin(deg2rad(α)))
    inputs.rho[1] = 1.225
end

# # 3. Integration and Results
# We use `simulate` to integrate the continuous-time state-space system and automatically extract the forces.
res = simulate(model, tspan; input_func=wagner_inputs!, solver=Tsit5())
CL_UVLM = [-cs[3] for cs in res.outputs.Cs];

# # 4. Analytical Solution: Wagner's Function (Jones Approximation)
# Jones Approximation of Wagners function is given by $w(s) = 1 - 0.165*exp(-0.0455*s) - 0.335*exp(-0.3*s)$
wagner(s) = 1.0 - 0.165*exp(-0.0455*s) - 0.335*exp(-0.3*s)
time_chord = res.t ./ (c/V) # Non-dimensional time (t*V/c)
s = 2.0 .* res.t ./ (c/V)       # Wagner's s (2*V*t/c)
CL_steady = 2 * π * deg2rad(α) # Thin airfoil theory steady lift
CL_analytical = CL_steady .* wagner.(s);

# # 5. Plotting
fig = Figure()
ax = Axis(fig[1, 1],
    title = "Sudden Acceleration (Wagner's Problem)",
    xlabel = L"t V / c",
    ylabel = L"C_L/C_{L,steady}")

lines!(ax, time_chord, CL_UVLM/CL_steady,
       label="UVLM (AeroPanels.jl)", linewidth=2)
lines!(ax, time_chord, CL_analytical/CL_steady,
       label="Wagner (Jones approx.)", linestyle=:dash, linewidth=2)
axislegend(ax, position = :rb)
fig

# # 6. Final check
println("Final CL (UVLM): ", CL_UVLM[end])
println("Final CL (Steady Theory): ", CL_steady)
println("Error: ", abs(CL_UVLM[end] - CL_steady)/CL_steady * 100, "%")
