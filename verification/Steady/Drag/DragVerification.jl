using AeroPanels
using StaticArrays
using ForwardDiff
using LaTeXStrings
using JSON

# --- Configuration ---
save_data = false  # Set to true to overwrite golden data, false to verify against it
plot_data = false
data_file = joinpath(@__DIR__, "drag_data.json")
plot_file = joinpath(@__DIR__, "AeroPanels.png")
tol = 1e-4


function calculate_aero_data(AR; nc=5, ns=60)
    chord = 1.0
    span = AR * chord / 2
    
    surf1 = AeroSurface2D(nc, ns, chord=(chord, chord), span=span, sweep=0.0)
    surf2 = Mirror(surf1, 2)
    
    props = AeroModelProperties(c=chord, b=2*span, S=2*span*chord)
    model = AeroModel2D([surf1, surf2], props)
    
    sol = AeroSolve(50.0, 2.0, model)
    cl = -sol.coefficientStability[3]
    cd = -sol.coefficientStability[1]
    k = cd / cl^2
    
    f(α) = -AeroSolve(50.0, α, model).coefficientStability[3]
    cla = ForwardDiff.derivative(f, 0.0) * 57.3
    
    return cla, k
end

# Range
ars = collect(4.0:2.0:40.0)
println("Calculating Drag data...")
results = [calculate_aero_data(ar) for ar in ars]
current_cla = [r[1] for r in results]
current_k = [r[2] for r in results]

if save_data
    data_to_persist = Dict(
        "ars" => ars,
        "cla" => current_cla,
        "k" => current_k
    )
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
    
    if !isapprox(current_cla, golden_data["cla"], atol=tol)
        error("Regression detected in CLa results!")
    end
    if !isapprox(current_k, golden_data["k"], atol=tol)
        error("Regression detected in K results!")
    end
    println("Verification successful! No regressions found.")
end

if plot_data
        
    using CairoMakie

    # Plotting settings
    set_theme!(fontsize = 24, font = "Latin Modern Math")

    fig = Figure(size=(1300, 600))

    # CLa Plot
    ax1 = Axis(fig[1, 1], title=L"Lift Curve Slope ($C_{L\alpha}$)",
            xlabel=L"AR", ylabel=L"C_{L\alpha} \text{ [1/rad]}")
    lines!(ax1, ars, current_cla, color=:blue, linewidth=3)
    scatter!(ax1, ars, current_cla, color=:blue, markersize=12)
    hlines!(ax1, [2π], color=:red, linestyle=:dash, label=L"2\pi", linewidth=2)
    ylims!(ax1, 3.5, 6.5)
    xlims!(ax1, 0, 40)

    # K Plot
    ax2 = Axis(fig[1, 2], title=L"Induced Drag Factor ($C_D/C_L^2$)",
            xlabel=L"AR", ylabel=L"C_D/C_L^2")
    lines!(ax2, ars, current_k, color=:black, linewidth=3)
    scatter!(ax2, ars, current_k, color=:black, markersize=12)
    lines!(ax2, ars, 1 ./ (π .* ars), color=:red, linestyle=:dash, label=L"1/(\pi AR)", linewidth=2)
    ylims!(ax2, 0, 0.08)
    xlims!(ax2, 0, 40)

    axislegend(ax1, position=:rb)
    axislegend(ax2, position=:rt)

    save(plot_file, fig)
    println("Plot updated at $plot_file")
end