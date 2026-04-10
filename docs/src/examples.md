# Examples

## Single-Firm ESRI with Components

```@doctest
using ESRI, SparseArrays, LinearAlgebra, Random
Random.seed!(42)

N = 30
W = sprand(N, N, 0.1) + 0.1I
info = IndustryInfo(rand(1:3, N), [true, false, true])
econ = ESRIEconomy(W, info)

res = esri(econ, 5; details = true, maxiter = 30, tol = 1e-3)
typeof(res), length(res.upstream), length(res.downstream)
# output
(ESRIResult{Float64}, 30, 30)
```

## Economy-Wide ESRI

```@doctest
scores = esri(econ; maxiter = 30, tol = 1e-3, threads = false)
length(scores) == 30 && all(0 .<= scores .<= 1)
# output
true
```

## Custom Final Weights and Combine Mode

```@doctest
weights = rand(Float32, N)
scores_up = esri(econ; final_weights = weights, combine = :upstream, maxiter = 25, tol = 1e-3)
length(scores_up), eltype(scores_up)
# output
(30, Float64)
```

## Single-Firm Partial Components

```@doctest
up_only = esri(econ, 3; components = :upstream, maxiter = 25, tol = 1e-3)
down_only = esri(econ, 3; components = :downstream, maxiter = 25, tol = 1e-3)
haskey(up_only, :upstream), haskey(down_only, :downstream)
# output
(true, true)
```

## Direct Shock Scenarios

```@doctest
psi = ones(N)
psi[1:3] .= 0.4
scenario = esri_shock(econ, psi; details = true, maxiter = 25, tol = 1e-3)
scenario.esri >= 0
# output
true
```

## Subset Computation

```@doctest
subset_scores = esri(econ; firm_indices = [2, 6, 9], maxiter = 20, tol = 1e-3, threads = false)
subset_scores[2] >= 0 && subset_scores[6] >= 0 && subset_scores[9] >= 0
# output
true
```

## Matrix-First Convenience Wrappers

```@doctest
value = compute_esri(W, info, 7; maxiter = 20, tol = 1e-3)
scenario_value = compute_esri_shock(W, info, psi; maxiter = 20, tol = 1e-3)
value >= 0 && scenario_value >= 0
# output
true
```
