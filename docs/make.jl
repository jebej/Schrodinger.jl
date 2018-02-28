using Documenter, Schrodinger, PyPlot

makedocs(
    sitename = "Schrodinger.jl",
    authors = "Jérémy Béjanin.",
    modules = [Schrodinger],
    linkcheck = true,
    clean = false,
    format = :html,
    html_prettyurls = false,
    pages = Any[
        "Home" => "index.md",
        "Manual" => Any[
            "man/gettingstarted.md",
            "man/quantumobjects.md",
            "man/working.md",
            "man/dynamics.md",
        ],
        "Examples" => Any[
            "examples/DRAG.md",
        ],
        "API" => Any[
            "api/quobj.md",
            "api/states.md",
            "api/operators.md",
            "api/functions.md",
        ],
    ]
)

deploydocs(
    repo = "github.com/jebej/Schrodinger.jl.git",
    julia = "0.6",
    target = "build",
    deps = nothing,
    make = nothing,
)
