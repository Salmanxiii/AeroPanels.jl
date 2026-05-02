using AeroPanels
using StaticArrays
using ForwardDiff
using JSON
using LaTeXStrings
using CairoMakie

# Set global theme
set_theme!(fontsize = 24, font = "DejaVu Sans")

# --- 1. Drag Verification Plot ---

function calculate_drag_data(AR; nc=5, ns=60)
    chord = 1.0
    span = AR * chord / 2
    surf1 = AeroSurface2D(nc, ns, chord=(chord, chord), span=span, sweep=0.0)
    surf2 = Mirror(surf1, 2)
    props = AeroModelProperties(c=chord, b=2*span, S=2*span*chord)
    model = AeroModel2D([surf1, surf2], props)
    f(α) = -AeroSolve(50.0, α, model).coefficientStability[3]
    cla = ForwardDiff.derivative(f, 0.0) * 57.3
    sol = AeroSolve(50.0, 2.0, model)
    cl = -sol.coefficientStability[3]
    cd = -sol.coefficientStability[1]
    k = cd / cl^2
    return cla, k
end

ars_drag = collect(4.0:4.0:40.0)
res_drag = [calculate_drag_data(ar) for ar in ars_drag]
cla_drag = [r[1] for r in res_drag]
k_drag = [r[2] for r in res_drag]

fig_drag = Figure(size=(1300, 600))
ax1 = Axis(fig_drag[1, 1], title=L"Lift Curve Slope ($C_{L\alpha}$)", xlabel=L"AR", ylabel=L"C_{L\alpha} \text{ [1/rad]}")
lines!(ax1, ars_drag, cla_drag, color=:blue, linewidth=3)
scatter!(ax1, ars_drag, cla_drag, color=:blue, markersize=12)
hlines!(ax1, [2π], color=:red, linestyle=:dash, label=L"2\pi", linewidth=2)
ylims!(ax1, 3.5, 6.5)
xlims!(ax1, 0, 40)

ax2 = Axis(fig_drag[1, 2], title=L"Induced Drag Factor ($C_D/C_L^2$)", xlabel=L"AR", ylabel=L"C_D/C_L^2")
lines!(ax2, ars_drag, k_drag, color=:black, linewidth=3)
scatter!(ax2, ars_drag, k_drag, color=:black, markersize=12)
lines!(ax2, ars_drag, 1 ./ (π .* ars_drag), color=:red, linestyle=:dash, label=L"1/(\pi AR)", linewidth=2)
ylims!(ax2, 0, 0.08)
xlims!(ax2, 0, 40)

save("docs/src/assets/AeroPanels.png", fig_drag)

# --- 2. Sweep Verification Plot ---

function calculate_sweep_cla(AR, sweep_deg; nc=5, ns=40)
    chord = 1.0
    span = AR * chord / 2
    surf1 = AeroSurface2D(nc, ns, chord=(chord, chord), span=span, sweep=deg2rad(sweep_deg))
    surf2 = Mirror(surf1, 2)
    props = AeroModelProperties(c=chord, b=2*span, S=2*span*chord)
    model = AeroModel2D([surf1, surf2], props)
    f(α) = -AeroSolve(50.0, α, model).coefficientStability[3]
    return ForwardDiff.derivative(f, 0.0) * 57.3
end

ars_sweep = collect(1.0:1.0:7.0)
sweeps = [0, 45, 60]
colors = [:blue, :red, :black]

fig_sweep = Figure(size=(800, 600))
ax_sweep = Axis(fig_sweep[1, 1], title="Sweep and Aspect Ratio Effect on CLa", xlabel=L"Aspect Ratio $AR$", ylabel=L"$C_{L\alpha}$")
xlims!(ax_sweep, 0, 8.5); ylims!(ax_sweep, 0, 7)
ax_sweep.xticks = 0:9; ax_sweep.yticks = 0:7

for (i, s) in enumerate(sweeps)
    res = [calculate_sweep_cla(ar, s) for ar in ars_sweep]
    lines!(ax_sweep, ars_sweep, res, label="$(s)° Sweep", color=colors[i], linewidth=2)
    scatter!(ax_sweep, ars_sweep, res, color=colors[i])
    val_inf = calculate_sweep_cla(100.0, s)
    lines!(ax_sweep, [8.9, 9.0], [val_inf, val_inf], color=colors[i], linewidth=2)
    text!(ax_sweep, 9.05, val_inf, text=string(round(val_inf, digits=2)), color=colors[i], align=(:left, :center), fontsize=14)
end
text!(ax_sweep, 9.05, 6.8, text="AR=100", color=:black, font=:bold, align=(:left, :center))
axislegend(ax_sweep, position=:lt)

save("docs/src/assets/Sweep_AeroPanels.png", fig_sweep)
