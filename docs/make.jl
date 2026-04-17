using Documenter
using ESRIcascade

DocMeta.setdocmeta!(ESRIcascade, :DocTestSetup, :(using ESRIcascade, SparseArrays, LinearAlgebra, Random); recursive = true)

makedocs(
    format = Documenter.HTML(edit_link = nothing),
    modules = [ESRIcascade],
    sitename = "ESRIcascade.jl",
    checkdocs = :exports,
    doctest = true,
    pages = [
        "Home" => "index.md",
        "Mathematical Model" => "theory.md",
        "Examples" => "examples.md",
        "Performance" => "performance.md",
        "Troubleshooting" => "troubleshooting.md",
        "API Reference" => "api.md",
    ],
)

if haskey(ENV, "DOCUMENTER_KEY")
    deploydocs(
        repo = "github.com/Devetak/ESRIcascade.jl.git",
        devbranch = "main",
    )
end
