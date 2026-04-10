# API Reference

## Core types

`IndustryInfo(industry_of_firm, essential_industry)` stores the industry index of each firm and one Boolean essentiality flag for each industry. The package assumes one-based industry ids. The constructor checks that each firm industry lies in `1:length(essential_industry)`.

`ESRIEconomy(W, info)` precomputes the upstream operator, the essential downstream operator, the non-essential downstream operator, row sums, column sums, total output, and firm count. Use this object when the network does not change across runs.

`ESRIResult` stores a scalar `esri` value and two vectors, `upstream` and `downstream`. The vectors are the converged health levels under one scenario.

## Main functions

`esri(econ; maxiter=100, tol=1e-2, verbose=false, threads=false, firm_indices=nothing, final_weights=nothing, combine=:min)` computes the default single-firm shock for all firms or for the firms selected by `firm_indices`. The default shock for firm `i` is the vector `\psi` with `\psi_i = 0` and `\psi_j = 1` for `j \neq i`.

`esri(econ, firm_idx; maxiter=100, tol=1e-2, verbose=false, details=false, components=:none, final_weights=nothing, combine=:min, shock=nothing)` computes one scenario and returns one scalar or one structured result.

`esri_shock(econ, shock; maxiter=100, tol=1e-2, verbose=false, details=false, components=:none, final_weights=nothing, combine=:min)` computes one scenario from an explicit capacity cap vector.

`compute_esri(W, info; ...)`, `compute_esri(W, info, firm_idx; ...)`, and `compute_esri_shock(W, info, shock; ...)` rebuild `ESRIEconomy` on each call and then dispatch to the corresponding `esri` method.

## Keyword semantics

`maxiter` is the maximum number of joint upstream and downstream iterations for one scenario.

`tol` is the stopping threshold for the infinity norm of the change between successive iterates.

`verbose=true` prints iteration information in single-scenario calls. In threaded economy-wide runs it is ignored and the package emits a warning.

`threads=true` applies only to economy-wide `esri(econ; ...)` and `compute_esri(W, info; ...)`. It parallelizes across shocked firms. It does not parallelize one fixed-point solve for one custom scenario.

`firm_indices` restricts the set of default firm shocks computed by `esri(econ; ...)`. Entries outside that set remain zero in the returned vector.

`details=true` forces a full `ESRIResult` return with both converged health vectors.

`components=:upstream` returns `(esri = value, upstream = vector)`. `components=:downstream` returns `(esri = value, downstream = vector)`. `components=:both` returns `ESRIResult`. `components=:none` returns only the scalar value.

`combine` selects the channel used in the final loss aggregation. `:min` uses `\min(u_i, d_i)`. `:upstream` uses `u_i`. `:downstream` uses `d_i`.

`final_weights` replaces the default output weights in the numerator of the ESRI reduction. The denominator stays equal to total output `\sum_i r_i`. It does not renormalize by `\sum_i w_i`.

`shock` is a firm-level capacity cap vector `\psi \in [0,1]^N`. `0` means full closure. `1` means no exogenous cap. Intermediate values impose partial capacity caps. In `esri(econ, firm_idx; shock=psi)`, the supplied `psi` is the whole scenario. The package does not apply an extra single-firm closure on top of it.

## Return values

`esri(econ; ...)` and `compute_esri(W, info; ...)` return a vector with one ESRI value per firm.

`esri(econ, firm_idx; ...)`, `esri_shock(econ, shock; ...)`, `compute_esri(W, info, firm_idx; ...)`, and `compute_esri_shock(W, info, shock; ...)` return one of three shapes. The scalar case returns one number. The partial-components case returns a named tuple with `esri` and the requested vector. The full-details case returns `ESRIResult`.

## Validation

The package throws `BoundsError` for an out-of-range firm index or an out-of-range entry in `firm_indices`. It throws `ArgumentError` for invalid `combine`, invalid `components`, duplicate `firm_indices`, empty `essential_industry`, or invalid industry ids. It throws `DimensionMismatch` if the matrix size does not match the metadata size, if `final_weights` has the wrong length, or if `shock` has the wrong length. It throws `DomainError` if `shock` contains values outside `[0,1]`, if `final_weights` contains negative or non-finite values, or if `W` contains negative or non-finite entries.

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
