# ESRI.jl

`ESRI.jl` computes firm-level economic systemic risk indicators from an input-output network.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/Devetak/ESRI.jl")
```

## Quick Start

```@doctest
using ESRI, SparseArrays, LinearAlgebra

N = 50
W = sprand(N, N, 0.08) + 0.1I
industry_ids = rand(1:4, N)
essential_industry = [true, true, false, false]
info = IndustryInfo(industry_ids, essential_industry)

econ = ESRIEconomy(W, info)
scores = esri(econ; maxiter = 40, tol = 1e-3)
length(scores), minimum(scores) >= 0, maximum(scores) <= 1
# output
(50, true, true)
```

## Documentation Map

- [Mathematical Model](@ref): equations and operator definitions used by the implementation.
- [Examples](@ref): doctested usage for single-firm, economy-wide, and custom shock analysis.
- [Performance Notes](@ref): sparse/dense behavior and threading guidance.
- [Troubleshooting](@ref): common errors and exact fixes.
- [API Reference](@ref): exported types and functions.

## Public API (Exported)

- `IndustryInfo`
- `ESRIEconomy`
- `esri`
- `esri_shock`
- `compute_esri`
- `compute_esri_shock`

For equations and algorithm details, see [Mathematical Model](@ref).
