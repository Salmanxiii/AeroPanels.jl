using Documenter
using AeroPanels
using Literate
using OrdinaryDiffEqTsit5
using StaticArrays
using GLMakie

# 1. Generate ALL Plots (Verification + Example Images)
# This script runs in the Main process where the environment is stable.
include("generate_verification_images.jl")
include("generate_example_images.jl")

# 2. Process Examples with Literate
example_dir = joinpath(@__DIR__, "..", "examples")
output_dir  = joinpath(@__DIR__, "src", "generated")

# Ensure output directory exists
mkpath(output_dir)

for file in readdir(example_dir)
    if endswith(file, ".jl")
        Literate.markdown(joinpath(example_dir, file), output_dir; 
            documenter=false, execute=true)
    end
end

# 3. Build Docs
makedocs(
    sitename = "AeroPanels.jl",
    modules  = [AeroPanels],
    doctest  = false,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    pages = [
        "Home" => "index.md",
        "Steady Theory" => "steady_theory.md",
        "Unsteady Theory" => "unsteady_theory.md",
        "Examples" => [
            "Steady Wing" => "generated/01_steady_wing.md",
            "Wagner Problem" => "generated/02_wagner_problem.md",
        ],
        "API Reference" => "api.md",
    ]
)

deploydocs(
    repo = "github.com/Salmanxiii/AeroPanels.jl.git",
    devbranch = "main",
)
