using Documenter, CollectiveSpins

pages = [
    "index.md",
    "installation.md",
    "geometry.md",
    "interaction.md",
    "effective_interaction.md",
    "system.md",
    "collective_modes.md",
    "reducedspin.md",
    "descriptions.md",
    "api.md"
]

makedocs(
    modules = [CollectiveSpins],
    checkdocs = :exports,
    format = Documenter.HTML(
        edit_link = nothing,
        canonical = "https://github.io/qojulia/CollectiveSpins.jl/",
        assets = [ asset("assets/favicon.png", class=:ico, islocal = true) ]),
    sitename = "CollectiveSpins.jl",
    pages = pages
    )

deploydocs(
    repo = "github.com/qojulia/CollectiveSpins.jl.git",
    )
