# API Reference

## Types

`IndustryInfo(industry_of_firm, essential_industry)` stores one industry id per firm and one Boolean essentiality flag per industry. Industry ids are one-based.

`ESRIEconomy(W, info)` caches the normalized upstream/downstream operators, row sums, column sums, total output, and firm count.

`ESRIResult` stores `esri`, `upstream`, and `downstream` for one scenario.

## Calls

`esri(econ; maxiter=100, tol=1e-2, verbose=false, threads=false, firm_indices=nothing, final_weights=nothing, combine=:min)` returns one ESRI value per firm. If `firm_indices` is set, only those default single-firm closures are computed and all other entries stay zero.

`esri(econ, firm_idx; maxiter=100, tol=1e-2, verbose=false, details=false, components=:none, final_weights=nothing, combine=:min, shock=nothing)` solves one scenario. By default it closes `firm_idx`. If `shock=psi` is given, `psi ∈ [0,1]^N` is the whole scenario.

`esri_shock(econ, shock; ...)` solves one scenario from an explicit capacity-cap vector `shock ∈ [0,1]^N`.

`compute_esri(...)` and `compute_esri_shock(...)` rebuild `ESRIEconomy` and dispatch to the matching `esri` method.

## Keywords

`maxiter` is the maximum number of upstream/downstream iterations for one scenario. `tol` is the infinity-norm stopping threshold.

`threads=true` applies only to economy-wide `esri(econ; ...)` and `compute_esri(W, info; ...)`. It parallelizes across shocked firms, not within one fixed-point solve.

`details=true` is shorthand for `components=:both`.

`components=:none` returns a scalar. `:upstream` returns `(esri=value, upstream=vector)`. `:downstream` returns `(esri=value, downstream=vector)`. `:both` returns `ESRIResult`.

`combine=:min` uses `min(upstream, downstream)` in the final reduction. `:upstream` and `:downstream` use that channel alone.

`final_weights` replaces the default output weights only in the numerator. The denominator stays `sum(row_sums)`.

`verbose=true` prints iteration logs for single-scenario calls. In threaded economy-wide runs it is ignored with a warning.

## Validation

`BoundsError` covers out-of-range firm ids. `ArgumentError` covers invalid `combine`, invalid `components`, duplicate `firm_indices`, empty `essential_industry`, and invalid industry ids. `DimensionMismatch` covers shape mismatches for `W`, `final_weights`, and `shock`. `DomainError` covers non-finite or negative `W`, non-finite or negative `final_weights`, and `shock` values outside `[0,1]`.

```@docs
IndustryInfo
ESRIEconomy
ESRIResult
```

```@docs
esri
esri_shock
compute_esri
compute_esri_shock
```
