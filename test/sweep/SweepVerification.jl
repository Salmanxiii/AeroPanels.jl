using AeroPanels
using StaticArrays
using ForwardDiff
using JSON
using Test

# --- Configuration ---
data_file = joinpath(@__DIR__, "sweep_data.json")
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

@testset "Sweep Verification" begin
    # Parameters
    ars = collect(1.0:1.0:7.0)
    sweeps = [0, 45, 60]

    current_results = Dict(string(s) => [calculate_cla(ar, s) for ar in ars] for s in sweeps)
    current_inf = Dict(string(s) => calculate_cla(100.0, s) for s in sweeps)

    if !isfile(data_file)
        error("Golden data file not found at $data_file.")
    end
    
    golden_data = JSON.parsefile(data_file)
    
    for s_str in keys(current_results)
        @test current_results[s_str] ≈ golden_data["results"][s_str] atol=tol
        @test current_inf[s_str] ≈ golden_data["inf_results"][s_str] atol=tol
    end
end
