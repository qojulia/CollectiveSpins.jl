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
    Documenter = Documenter.HTML(
        edit_link = nothing,
        canonical = "https://github.io/qojulia/CollectiveSpins.jl/",
        assets = [ asset("assets/favicon.png", class=:ico, islocal = true) ]),
    sitename = "CollectiveSpins.jl",
    pages = pages
    )