using Documenter, CollectiveSpins

pages = [
    "index.md",
    "installation.md",
    "geometry.md",
    "interaction.md",
    "effective_interaction.md",
    "system.md",
    "descriptions.md",
    "api.md"
]

makedocs(
    modules = [CollectiveSpins],
    checkdocs = :exports,
    format = :html,
    sitename = "CollectiveSpins.jl",
    pages = pages
    )

deploydocs(
    repo = "github.com/bastikr/CollectiveSpins.jl.git",
    target = "build",
    julia = "0.6",
    deps = nothing,
    make = nothing
)
