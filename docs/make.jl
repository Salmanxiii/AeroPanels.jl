using Documenter
using AeroPanels

makedocs(
    sitename = "AeroPanels.jl",
    modules  = [AeroPanels],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    pages = [
        "Home" => "index.md",
        "Steady Theory" => "steady_theory.md",
        "Unsteady Theory" => "unsteady_theory.md",
        "Examples" => "examples.md",
        "API Reference" => "api.md",
    ]
)

deploydocs(
    repo = "github.com/Salmanxiii/AeroPanels.jl.git",
    devbranch = "main",
)
