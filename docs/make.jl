using Documenter, Schrodinger

makedocs(
    sitename = "Schrodinger.jl",
    authors = "Jérémy Béjanin.",
    modules = [Schrodinger],
    linkcheck = true,
    clean = false,
    format = :html,
    html_prettyurls = true,
    pages = Any[
        "Home" => "index.md",
        "Manual" => Any[
            "man/gettingstarted.md",
            "man/quantumobjects.md",
            "man/working.md",
            "man/dynamics.md",
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
    julia = "0.5",
    target = "build",
    deps = nothing,
    make = nothing,
    )
