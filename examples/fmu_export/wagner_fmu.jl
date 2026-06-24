using AeroPanels
using StaticArrays
using JuliaC
using XML
using CodeTracking

# --- 1. Model Definition ---
# This function defines the aerodynamic model that we want to export as an FMU.
# It MUST return an UnsteadyAeroModel2D.
function create_wagner_model()
    c, V, α = 1.0, 1.0, 2.0
    # High-aspect ratio wing to approximate 2D behavior
    surf = AeroSurface2D(20, 1, chord=(c, c), span=100.0)
    props = AeroModelProperties(c=c, b=100.0, S=100.0)
    
    # We use higher fidelity for the FMU example to match the reference results.
    # Reference plot uses nWake=500, wakeLength=20.0. 
    # We use nWake=200, wakeLength=15.0 for a balance of speed and accuracy.
    model = UnsteadyAeroModel2D([surf], props, V, nWake=200, wakeLength=15.0)
    return model
end

fmu_name = "WagnerWingFMU"
fmu_dir = @__DIR__
fmu_path = joinpath(fmu_dir, fmu_name * ".fmu")

# --- 2. Build the FMU ---
# This step compiles the AeroPanels model and the FMI wrapper into a single .fmu file.
# It uses JuliaC for static compilation and a C-shim for Windows DLL loading.
println("Building FMU: $fmu_name...")
BuildFMU(create_wagner_model, fmu_dir, fmu_name=fmu_name, clean=true)

println("--- SUCCESS ---")
println("FMU created at: $fmu_path")
println("")
println("Note on Testing:")
println("The exported FMU uses FMI 2.0. You can test it using tools like:")
println("1. FMI.jl: loadFMU(\"$fmu_name.fmu\")")
println("2. Python fmpy: fmpy.simulate_fmu(\"$fmu_name.fmu\")")
println("")
println("Important: On Windows, loading a compiled Julia FMU into another Julia process")
println("can cause process-wide Julia state conflicts (signal 22). It is recommended")
println("to test the FMU in an environment without Julia (e.g., Python or Modelica)")
println("or in a separate, clean Julia session that doesn't have AeroPanels pre-loaded.")
