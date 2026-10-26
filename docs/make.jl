using Documenter
using GNSSSignals

makedocs(
    sitename = "GNSSSignals.jl",
    modules = [GNSSSignals],
    authors = "JuliaGNSS",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://JuliaGNSS.github.io/GNSSSignals.jl",
    ),
    pages = [
        "Home" => "index.md",
        "Usage" => "usage.md",
        "API Reference" => [
            "Signals" => "api.md",
            "Bands and Time Systems" => "bands_and_time_systems.md",
            "Codes" => "codes.md",
        ],
    ],
    checkdocs = :exports,
)

deploydocs(
    repo = "github.com/JuliaGNSS/GNSSSignals.jl.git",
    devbranch = "master",
    push_preview = true,
)
