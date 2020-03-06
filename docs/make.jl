using Documenter, Schrodinger, SparseArrays, PyPlot

DocMeta.setdocmeta!(Schrodinger, :DocTestSetup, :(using Schrodinger, LinearAlgebra, SparseArrays); recursive=true)

makedocs(
    sitename = "Schrodinger.jl",
    authors = "Jérémy Béjanin.",
    modules = [Schrodinger],
    format = Documenter.HTML(prettyurls = false),
    linkcheck = true,
    clean = false,
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "man/gettingstarted.md",
            "man/quantumobjects.md",
            "man/working.md",
            "man/dynamics.md",
        ],
        "Examples" => [
            "examples/driven_TLS.md",
            "examples/DRAG.md",
        ],
        "API" => [
            "api/quobj.md",
            "api/states.md",
            "api/operators.md",
            "api/functions.md",
        ],
    ]
)

deploydocs(
    repo = "github.com/jebej/Schrodinger.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
