using Documenter, CollectiveSpins

pages = [
    "index.md",
    "installation.md",
    "geometry.md",
    "interaction.md",
    "effective_interaction.md",
    "system.md",
    "descriptions.md"
]

makedocs(
    modules = [CollectiveSpins],
    checkdocs = :exports,
    format = :html,
    sitename = "CollectiveSpins.jl",
    pages = pages
    )
