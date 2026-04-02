# Troubleshooting

## `DimensionMismatch` from `ESRIEconomy`

Cause: `weight_matrix` is not square or does not match the number of firms in `IndustryInfo`.

Fix:

```julia
size(W, 1) == size(W, 2) == length(info.industry_of_firm)
```

## `ArgumentError` from `IndustryInfo`

Cause: `industry_of_firm` contains values outside `1:length(essential_industry)` or `essential_industry` is empty.

Fix:

```julia
minimum(industry_ids) >= 1 && maximum(industry_ids) <= length(essential_industry)
```

## Invalid `combine` or `components`

Cause: unsupported symbols.

Fix:

- `combine ∈ (:min, :upstream, :downstream)`
- `components ∈ (:none, :upstream, :downstream, :both)`

## Shock validation failures

Cause: `shock` has wrong length or entries outside `[0,1]`.

Fix:

```julia
length(shock) == econ.n
all(0 .<= shock .<= 1)
```

## Unexpected zeros in subset runs

When using `firm_indices`, ESRI is computed only for selected firms. Non-selected entries remain zero by design.

## Slow repeated calls

If you call `compute_esri` repeatedly with the same `W` and `info`, prebuild `econ = ESRIEconomy(W, info)` and call `esri(econ, ...)` / `esri_shock(econ, ...)` to avoid repeated preprocessing.
