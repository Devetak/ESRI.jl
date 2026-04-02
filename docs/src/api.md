# API Reference

## Constructors

- `IndustryInfo(industry_of_firm::AbstractVector{<:Integer}, essential_industry::AbstractVector{Bool})`
- `ESRIEconomy(weight_matrix, info::IndustryInfo)`

## Main Functions

- `esri(econ::ESRIEconomy; maxiter=100, tol=1e-2, verbose=false, threads=false, firm_indices=nothing, final_weights=nothing, combine=:min)`
- `esri(econ::ESRIEconomy, firm_idx::Integer; maxiter=100, tol=1e-2, verbose=false, details=false, components=:none, final_weights=nothing, combine=:min, shock=nothing)`
- `esri_shock(econ::ESRIEconomy, shock; maxiter=100, tol=1e-2, verbose=false, details=false, components=:none, final_weights=nothing, combine=:min)`
- `compute_esri(weight_matrix, info; kwargs...)`
- `compute_esri(weight_matrix, info, firm_idx::Integer; kwargs...)`
- `compute_esri_shock(weight_matrix, info, shock; kwargs...)`

## Return Shapes

- Economy-wide `esri(econ; ...)` / `compute_esri(W, info; ...)` return `Vector{T}` with one ESRI value per firm.
- Single-firm/scenario calls return:
  - scalar `T` when `components = :none` and `details = false`
  - named tuple for one component (`:upstream` or `:downstream`)
  - `ESRIResult{T}` when `details = true` or `components = :both`

## Validation and Errors

- `BoundsError`: invalid firm index / out-of-range `firm_indices`.
- `ArgumentError`: invalid `combine`, `components`, duplicate `firm_indices`, invalid industry ids.
- `DimensionMismatch`: size mismatch in matrix/info, `final_weights`, or `shock`.
- `DomainError`: shock values outside `[0, 1]`.

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
