# API Reference

## Types

`IndustryInfo(industry_of_firm)` infers the industry count from the largest one-based id and classifies every supplier-customer industry pair as essential.

`IndustryInfo(industry_of_firm, essential_industry)` stores one one-based industry id per firm and one Boolean essentiality flag per supplying industry. It is a backwards-compatible shorthand: an essential industry is classified as `2` for every customer industry, and a non-essential industry as `1` for every customer industry.

`IndustryInfo(industry_of_firm, input_classification)` accepts a square integer matrix whose rows are supplying industries and columns are customer industries. Entries use `2` for essential, `1` for non-essential, and `0` for no downstream production effect. Matrix indices and `industry_of_firm` must use the same one-based industry numbering. Input links classified as `0` still remain available to upstream propagation. Classification values are trusted rather than scanned during construction.

`ihs_input_classification()` returns the bundled 616 by 616 IHS matrix and `ihs_industry_codes()` returns its label order. They are opt-in defaults; `IndustryInfo(industry_of_firm)` remains the all-essential compatibility default.

`ESRIEconomy(W, info)` caches the normalized upstream/downstream operators, row sums, column sums, total output, and firm count.

`ESRIResult` stores `esri`, `upstream`, and `downstream` for one scenario.

## Calls

`esri(econ; maxiter=100, tol=1e-2, verbose=false, threads=false, firm_indices=nothing, final_weights=nothing, combine=:min)` returns one ESRI value per firm. If `firm_indices` is set, only those default single-firm closures are computed and all other entries stay zero.

`esri(econ, firm_idx; maxiter=100, tol=1e-2, verbose=false, details=false, components=:none, final_weights=nothing, combine=:min, shock=nothing)` solves one scenario. By default it closes `firm_idx`. If `shock=psi` is given, then `psi ∈ [0,1]^N` is the whole scenario: `psi[i] = 1` means no exogenous shock to firm `i`, `psi[i] = 0` means firm `i` is closed, and intermediate values cap firm `i` at partial capacity.

`esri_shock(econ, shock; ...)` solves one scenario from an explicit capacity-cap vector `shock ∈ [0,1]^N`.

`compute_esri(...)` and `compute_esri_shock(...)` rebuild `ESRIEconomy` and dispatch to the matching `esri` method.

## Keywords

`maxiter` is the maximum number of upstream/downstream iterations for one scenario. `tol` is the infinity-norm stopping threshold.

`threads=true` applies only to economy-wide `esri(econ; ...)` and `compute_esri(W, info; ...)`. It parallelizes across shocked firms, not within one fixed-point solve.

`details=true` is shorthand for `components=:both`.

`components=:none` returns a scalar. `:upstream` returns `(esri=value, upstream=vector)`. `:downstream` returns `(esri=value, downstream=vector)`. `:both` returns `ESRIResult`.

`combine=:min` uses `min(upstream, downstream)` in the final reduction. `:upstream` and `:downstream` use that channel alone.

`final_weights` replaces the default output weights in the final reduction. With the default weights, the scalar is a share of total output. With custom weights, the scalar is the weighted loss divided by `sum(final_weights)`.

`verbose=true` prints iteration logs for single-scenario calls. In threaded economy-wide runs it is ignored with a warning.

## Validation

`BoundsError` covers out-of-range firm ids. `ArgumentError` covers invalid `combine`, invalid `components`, duplicate `firm_indices`, empty classifications, and invalid industry ids. `DimensionMismatch` covers shape mismatches for `W`, `input_classification`, `final_weights`, and `shock`. `DomainError` covers non-finite or negative `W`, non-finite or negative `final_weights`, and `shock` values outside `[0,1]`.

```@docs
IndustryInfo
ESRIEconomy
ESRIResult
ihs_input_classification
ihs_industry_codes
```

```@docs
esri
esri_shock
compute_esri
compute_esri_shock
```
