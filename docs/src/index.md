# ESRI.jl

`ESRI.jl` computes firm-level ESRI values on a directed firm-to-firm supply network.

The package takes a square weight matrix `W`, a firm-to-industry map, and an industry-level essentiality vector. It builds a reusable `ESRIEconomy` object and then solves either one shock per firm or one custom shock vector.

The package is built around three exported entry points. `ESRIEconomy(W, info)` precomputes the linear operators and baseline weights. `esri(econ; ...)` computes one default firm shock for each selected firm. `esri_shock(econ, psi; ...)` computes one custom scenario with capacity cap vector `psi`.

## Installation

```julia
using Pkg
Pkg.add("ESRI")
```

For local documentation builds, run

```julia
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

## Quick start

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

The main performance rule is simple. Build `ESRIEconomy` once. Reuse it for all runs on the same network.

The theory page states the operators and fixed-point equations used by the package. The examples page states the calling patterns. The API page states the exact keyword behavior and return shapes.
