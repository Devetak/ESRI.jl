# ESRI.jl

`ESRI.jl` is a package for computing the Economic Systemic Risk Index for firms in an economy.

Inputs are a square weight matrix `W` representing the supply chain, one industry id per firm, and one Boolean essentiality flag per industry.

`ESRIEconomy(W, info)` caches the operators and weights for further reuse. `esri(econ; ...)` computes the default single-firm closures. `esri_shock(econ, psi; ...)` computes one explicit scenario from a firm-level capacity-cap vector `psi`, where `1` means unaffected and `0` means closed.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/Devetak/ESRI.jl")
```

## Quick start

```@doctest
using ESRI, SparseArrays, LinearAlgebra

N = 50
W = sprand(N, N, 0.08) + 0.1I
info = IndustryInfo(rand(1:4, N), [true, true, false, false]) # industry 1 and 2 are essential

econ = ESRIEconomy(W, info) # set up the economy
scores = esri(econ; maxiter = 40, tol = 1e-3) # compute ESRI for each firm
length(scores), minimum(scores) >= 0, maximum(scores) <= 1
# output
(50, true, true)
```

Build `ESRIEconomy` once and reuse it on the same network.
