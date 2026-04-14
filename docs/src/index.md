# ESRI.jl

`ESRI.jl` is a package for computing the Economic Systemic Risk Index for firms in an economy.

Inputs are a square weight matrix `W` representing the supply chain, one industry id per firm, and one Boolean essentiality flag per industry.

`ESRIEconomy(W, info)` caches the operators and weights for further reuse. `esri(econ; ...)` computes the default single-firm closures. `esri_shock(econ, psi; ...)` computes one explicit scenario from a firm-level capacity-cap vector `psi`, where `1` means unaffected and `0` means closed.

In simple terms, `psi[i]` says how much of its normal capacity firm `i` is allowed to use in the scenario you want to study. Use it to describe shocks like a plant shutdown, an energy shortage, a port disruption, sanctions, or a sector-wide restriction. ESRI then shows how that local shock can spread through suppliers and customers across the wider economy.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/Devetak/ESRI.jl")
```

## Quick start

```@doctest
using ESRI, SparseArrays

N = 1_000
W = sprand(N, N, 0.01)
W[1:N+1:end] .= 0
info = IndustryInfo(rand(1:4, N), [true, true, false, false]) # industry 1 and 2 are essential

econ = ESRIEconomy(W, info) # set up the economy
scores = esri(econ; maxiter = 40, tol = 1e-3) # compute ESRI for each firm
nothing
```

![Histogram of example ESRI scores](assets/scores_hist.svg)

The histogram is usually the first thing to look at. If most firms sit near zero, most single-firm failures have limited economy-wide spillovers. A longer right tail means some firms create much larger losses when they fail.

Build `ESRIEconomy` once and reuse it on the same network.
