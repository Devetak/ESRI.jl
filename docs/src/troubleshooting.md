# Troubleshooting

If `ESRIEconomy(W, info)` throws `DimensionMismatch`, then `W` is not square or its size does not match the number of firms in `info`.

If `IndustryInfo(industry_of_firm, essential_industry)` throws `ArgumentError`, then `essential_industry` is empty or some firm industry id lies outside `1:length(essential_industry)`.

If `esri(econ, firm_idx; ...)` throws `BoundsError`, then `firm_idx` is out of range. If `esri(econ; firm_indices=...)` throws `BoundsError`, then one entry of `firm_indices` is out of range.

If `combine` is invalid, the package accepts only `:min`, `:upstream`, and `:downstream`.

If `components` is invalid, the package accepts only `:none`, `:upstream`, `:downstream`, and `:both`.

If `shock` throws `DimensionMismatch`, then its length is not `econ.n`. If `shock` throws `DomainError`, then at least one value is not finite or does not lie in `[0,1]`.

If `final_weights` throws `DimensionMismatch`, then its length is not `econ.n`. If it throws `DomainError`, then at least one value is negative or not finite.

If an economy-wide run with `firm_indices` returns zeros outside the selected set, this is the intended result shape.

If `esri(econ, firm_idx; shock=psi)` does not match the default single-firm shock for `firm_idx`, this is expected. Once `shock=psi` is supplied, `psi` is the scenario.
