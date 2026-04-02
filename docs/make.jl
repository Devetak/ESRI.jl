using Documenter
using ESRI

DocMeta.setdocmeta!(ESRI, :DocTestSetup, :(using ESRI, SparseArrays, LinearAlgebra, Random); recursive = true)

makedocs(
    modules = [ESRI],
    sitename = "ESRI.jl",
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
        repo = "github.com/Devetak/ESRI.jl.git",
        devbranch = "main",
    )
end
