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
length(scores) == N && all(0 .<= scores .<= 1)
# output
true
```

## Single firm with full details

```@doctest
res = esri(econ, 5; details = true, maxiter = 30, tol = 1e-3)
typeof(res), length(res.upstream), length(res.downstream)
# output
(ESRIResult{Float64}, 30, 30)
```

## Requested components only

```@doctest
up_only = esri(econ, 3; components = :upstream, maxiter = 25, tol = 1e-3)
down_only = esri(econ, 3; components = :downstream, maxiter = 25, tol = 1e-3)
haskey(up_only, :upstream), haskey(down_only, :downstream)
# output
(true, true)
```

## Combine modes

```@doctest
score_min = esri(econ, 4; combine = :min, maxiter = 25, tol = 1e-3)
score_up = esri(econ, 4; combine = :upstream, maxiter = 25, tol = 1e-3)
score_down = esri(econ, 4; combine = :downstream, maxiter = 25, tol = 1e-3)
all(>=(0), [score_min, score_up, score_down])
# output
true
```

## Custom final weights

```@doctest
weights = rand(Float32, N)
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

The package takes a capacity cap vector `psi`. The entry `0.0` is a full closure. The entry `0.5` is a half-closure. The entry `1.0` leaves the firm uncapped.

```@doctest
psi = ones(N)
psi[1] = 0.0
psi[2] = 0.0
psi[3] = 0.5
psi[4] = 0.5
psi[5] = 0.8

scenario = esri_shock(econ, psi; details = true, maxiter = 25, tol = 1e-3)
scenario.esri >= 0 && all(scenario.upstream .<= psi .+ 1e-8) && all(scenario.downstream .<= psi .+ 1e-8)
# output
true
```

## Single-firm call with explicit shock

When `shock=psi` is supplied, the scenario is `psi`. The package does not add a second closure on `firm_idx`.

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
scenario_value = compute_esri_shock(W, info, psi; maxiter = 20, tol = 1e-3)
value >= 0 && scenario_value >= 0
# output
true
```
