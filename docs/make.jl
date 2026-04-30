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
        "Theory" => "theory.md",
        "Examples" => "examples.md",
        "API Reference" => "api.md",
    ]
)

deploydocs(
    repo = "github.com/Julius/AeroPanels.jl.git",
    devbranch = "main",
)
