using AeroPanels
using StaticArrays
using ForwardDiff
using JSON
using Test

# --- Configuration ---
data_file = joinpath(@__DIR__, "drag_data.json")
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

@testset "Drag Verification" begin
    # Range
    ars = collect(4.0:2.0:40.0)
    results = [calculate_aero_data(ar) for ar in ars]
    current_cla = [r[1] for r in results]
    current_k = [r[2] for r in results]

    if !isfile(data_file)
        error("Golden data file not found at $data_file.")
    end
    
    golden_data = JSON.parsefile(data_file)
    
    @test current_cla ≈ golden_data["cla"] atol=tol
    @test current_k ≈ golden_data["k"] atol=tol
end
