# API Reference

## Types

`IndustryInfo(industry_of_firm)` infers the industry count from the largest
one-based id. It classifies every supplier-customer pair as a linear input
(`1`).

`IndustryInfo(industry_of_firm, essential_industry)` stores one one-based
industry id per firm and one Boolean essentiality flag per supplying industry.
It is a backwards-compatible shorthand. An essential industry receives class
`2` for every customer industry. A linear industry receives class `1` for
every customer industry.

`IndustryInfo(industry_of_firm, input_classification)` accepts a square
integer-code matrix. `Matrix{UInt8}` is the standard form. A Boolean matrix
raises `ArgumentError` because its meaning is ambiguous. Rows represent
supplying industries. Columns represent customer industries. The codes are `2`
for essential input and `1` for linear input. Code `0` keeps the input in `W`
and the upstream operator and contributes to the linear-input denominator.

Matrix indices and `industry_of_firm` use the same one-based numbering. Class-0
links remain in `W` and the upstream operator. Class-1 links use the linear
downstream operator. Class-2 links use the essential downstream operator.
Supply an integer-code matrix such as `Matrix{UInt8}` to represent all three
states. The constructor trusts classification values during construction.
Validate external data before passing it in. The legacy `essential_industry`
and `essential_firm` fields remain meaningful for the Boolean-vector form. For
a customer-specific matrix, they are true when the supplier-industry row
consists entirely of `2`.

`ihs_input_classification()` returns the bundled 616 × 616 IHS matrix.
`ihs_industry_codes()` returns its label order. The bundled classification is
opt-in. `IndustryInfo(industry_of_firm)` uses the purely linear baseline. The
bundled IHS classification uses CC BY 4.0 and derives from the
[Pichler et al. version 3 archive](https://doi.org/10.5281/zenodo.7113762).

`ESRIEconomy(W, info)` caches the normalized upstream and downstream operators,
row sums, column sums, total output, and firm count. `W[i, j]` is the
nonnegative supply from supplier `i` to customer `j`.

`ESRIResult` stores `esri`, `upstream`, and `downstream` for one scenario.

## Calls

`esri(econ; ...)` returns one ESRI value per firm. Set `firm_indices` to select
default single-firm closures. The remaining entries stay zero.

`esri(econ, firm_idx; ...)` solves one scenario. The default scenario closes
`firm_idx`. Supply `shock=psi` to provide the complete scenario. Here
``\psi \in [0,1]^N``. ``\psi_i = 1`` makes full capacity available,
``\psi_i = 0`` closes firm ``i``, and intermediate values set a partial
capacity limit.

`esri_shock(econ, shock; ...)` solves one scenario from an explicit capacity-cap
vector in ``[0,1]^N``.

`compute_esri(W, info; ...)` rebuilds `ESRIEconomy` and dispatches to
economy-wide `esri(econ; ...)`. `compute_esri(W, info, firm_idx; ...)`
dispatches to the single-scenario overload with `details`, `components`,
`final_weights`, `combine`, and `shock`.

`compute_esri_shock(W, info, shock; ...)` rebuilds `ESRIEconomy` and dispatches
to `esri_shock`.

## Keywords

`maxiter` is the maximum number of upstream and downstream iterations for one
scenario. Use a positive value. `tol` is the infinity-norm stopping threshold.
The solver returns the final iterate when it reaches `maxiter`.

`threads=true` applies to economy-wide `esri(econ; ...)` and
`compute_esri(W, info; ...)`. It uses the available Julia threads across
shocked firms. With one Julia thread, the call runs serially. Each
single-scenario solve uses one thread.

See [Performance](performance.md) for the submitted runtime comparison.

`details=true` is shorthand for `components=:both`.

`components=:none` returns a scalar. `:upstream` returns
`(esri=value, upstream=vector)`. `:downstream` returns
`(esri=value, downstream=vector)`. `:both` returns `ESRIResult`.

Component vectors contain remaining-production fractions in `[0,1]`.
Upstream values reflect customer losses. Downstream values reflect input
shortages.

`combine=:min` uses `min(upstream, downstream)` in the final reduction.
`:upstream` and `:downstream` use that channel alone.

`final_weights` replaces the default output weights in the final reduction.
The default weights produce a share of total output. Custom weights produce a
weighted loss divided by `sum(final_weights)`. A zero total weight returns `0`.

`verbose=true` prints iteration logs for single-scenario calls. With multiple
Julia threads, `threads=true` and `verbose=true` emit a warning and disable the
progress display.

## Validation

- `BoundsError` identifies a firm id. Use ids in `1:econ.n`.
- `ArgumentError` identifies unsupported `combine` or `components` values,
  duplicate `firm_indices`, an empty classification, or an industry id outside
  `1:nindustries`.
- `DimensionMismatch` identifies shape errors in `W`, `input_classification`,
  `final_weights`, or `shock`.
- `DomainError` identifies value-domain errors. `W` and `final_weights` require
  finite values `>= 0`. `shock` requires values in `[0,1]`.

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
