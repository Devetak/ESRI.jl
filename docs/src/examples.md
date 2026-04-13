# Examples

## Economy-wide run

```@doctest
using ESRI, SparseArrays, LinearAlgebra, Random
Random.seed!(42)

N = 30
W = sprand(N, N, 0.1) + 0.1I
info = IndustryInfo(rand(1:3, N), [true, false, true])
econ = ESRIEconomy(W, info)

scores = esri(econ; maxiter = 30, tol = 1e-3, threads = false)
```

## Single firm with full details

```@doctest
res = esri(econ, 5; details = true, maxiter = 20, tol = 1e-3)
typeof(res), length(res.upstream), length(res.downstream)
# (ESRIResult{Float64}, 30, 30)
```

## Requested components only

```@doctest
up_only = esri(econ, 3; components = :upstream, maxiter = 25, tol = 1e-3)
down_only = esri(econ, 3; components = :downstream, maxiter = 25, tol = 1e-3)
```

## Custom final weights

```@doctest
weights = rand(Float32, N) # for example number of employees
scores_up = esri(econ; final_weights = weights, combine = :upstream, maxiter = 25, tol = 1e-3)
length(scores_up), eltype(scores_up)
# output
(30, Float64)
```

## Subset of default firm shocks

```@doctest
subset_scores = esri(econ; firm_indices = [2, 6, 9], maxiter = 20, tol = 1e-3, threads = false)
subset_scores[2] >= 0 && subset_scores[6] >= 0 && subset_scores[9] >= 0 && all(subset_scores[[1,3,4,5,7,8,10:end]] .== 0)
# output
true
```

## Custom shock vector

```@doctest
psi = ones(N)
psi[1] = 0.0
psi[2] = 0.0
psi[3] = 0.5
psi[4] = 0.5
psi[5] = 0.8

scenario = esri_shock(econ, psi; details = true, maxiter = 25, tol = 1e-3)
```

## Single-firm call with explicit shock

```@doctest
psi2 = ones(N)
psi2[1:3] .= 0.4

res1 = esri(econ, 7; shock = psi2, details = true, maxiter = 25, tol = 1e-3)
res2 = esri_shock(econ, psi2; details = true, maxiter = 25, tol = 1e-3)
res1.esri ≈ res2.esri
# output
true
```

## Matrix-first wrappers

```@doctest
value = compute_esri(W, info, 7; maxiter = 20, tol = 1e-3)
```
