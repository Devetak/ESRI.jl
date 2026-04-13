# ESRI.jl

`ESRI.jl` computes firm-level ESRI values on directed supply networks.

Inputs are a square weight matrix `W`, one industry id per firm, and one Boolean essentiality flag per industry. This package keeps the paper's upstream/downstream ESRI structure, but it narrows the inputs to that simpler contract and treats `psi` as a capacity-cap vector in `[0,1]^N`.

`ESRIEconomy(W, info)` caches the operators and weights. `esri(econ; ...)` computes the default single-firm closures. `esri_shock(econ, psi; ...)` computes one explicit scenario.

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
info = IndustryInfo(rand(1:4, N), [true, true, false, false])

econ = ESRIEconomy(W, info)
scores = esri(econ; maxiter = 40, tol = 1e-3)
length(scores), minimum(scores) >= 0, maximum(scores) <= 1
# output
(50, true, true)
```

Build `ESRIEconomy` once and reuse it on the same network. The theory page gives the equations, the API page gives keyword behavior and return shapes, and the examples page shows the main call patterns.
