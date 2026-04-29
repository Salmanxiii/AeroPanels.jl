using AeroPanels
using StaticArrays
using ForwardDiff
using JSON
using LaTeXStrings

# --- Configuration ---
save_data = false  # Set to true to overwrite golden data, false to verify against it
save_plot = false
data_file = joinpath(@__DIR__, "sweep_data.json")
plot_file = joinpath(@__DIR__, "AeroPanels.png")
tol = 1e-4


function calculate_cla(AR, sweep_deg; nc=5, ns=40)
    chord = 1.0
    span = AR * chord / 2
    sweep = deg2rad(sweep_deg)
    
    surf1 = AeroSurface2D(nc, ns, chord=(chord, chord), span=span, sweep=sweep)
    surf2 = Mirror(surf1, 2)
    
    props = AeroModelProperties(c=chord, b=2*span, S=2*span*chord)
    model = AeroModel2D([surf1, surf2], props)
    
    f(α) = -AeroSolve(50.0, α, model).coefficientStability[3]
    return ForwardDiff.derivative(f, 0.0) * 57.3
end

# Parameters
ars = collect(1.0:1.0:7.0)
sweeps = [0, 45, 60]
colors = [:blue, :red, :black]

println("Calculating Sweep data...")
current_results = Dict(string(s) => [calculate_cla(ar, s) for ar in ars] for s in sweeps)
current_inf = Dict(string(s) => calculate_cla(100.0, s) for s in sweeps)

data_to_persist = Dict(
    "results" => current_results,
    "inf_results" => current_inf,
    "ars" => ars
)

if save_data
    open(data_file, "w") do f
        JSON.print(f, data_to_persist, 4)
    end
    println("Golden data saved to $data_file")
else
    if !isfile(data_file)
        error("Golden data file not found at $data_file. Run with save_data = true first.")
    end
    
    golden_data = JSON.parsefile(data_file)
    println("Verifying against golden data...")
    
    for s_str in keys(current_results)
        # Verify main curves
        if !isapprox(current_results[s_str], golden_data["results"][s_str], atol=tol)
            error("Regression detected in results for sweep $(s_str)!")
        end
        # Verify inf results
        if !isapprox(current_inf[s_str], golden_data["inf_results"][s_str], atol=tol)
            error("Regression detected in infinite AR results for sweep $(s_str)!")
        end
    end
    println("Verification successful! No regressions found.")
end

# --- Plotting ---
if save_plot
    using CairoMakie
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1], title="Sweep and Aspect Ratio Effect on CLa",
            xlabel=L"Aspect Ratio $AR$", ylabel=L"$C_{L\alpha}$")

    xlims!(ax, 0, 8.5)
    ylims!(ax, 0, 7)
    ax.xticks = 0:9
    ax.yticks = 0:7

    for (i, sweep) in enumerate(sweeps)
        s_str = string(sweep)
        res = current_results[s_str]
        lines!(ax, ars, res, label="$(s_str)° Sweep", color=colors[i], linewidth=2)
        scatter!(ax, ars, res, color=colors[i])
        
        val_inf = current_inf[s_str]
        lines!(ax, [8.9, 9.0], [val_inf, val_inf], color=colors[i], linewidth=2)
        text!(ax, 9.05, val_inf, text=string(round(val_inf, digits=2)), 
            color=colors[i], align=(:left, :center), fontsize=14)
    end

    text!(ax, 9.05, 6.8, text="AR=100", color=:black, font=:bold, align=(:left, :center))
    axislegend(ax, position=:lt)
    save(plot_file, fig)
    println("Plot updated at $plot_file")
end