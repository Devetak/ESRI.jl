# Troubleshooting

- `ESRIEconomy(W, info)` raises `DimensionMismatch`: use a square `W` with one
  row and one column per firm in `info`.
- `IndustryInfo(industry_of_firm, essential_industry)` raises `ArgumentError`:
  provide one Boolean flag per industry and use ids from
  `1:length(essential_industry)`.
- `IndustryInfo(industry_of_firm, input_classification)`: use a square
  `Matrix{UInt8}` with at least one row and codes `0`, `1`, and `2`. Use matching
  one-based industry ids. Validate external matrix data before construction.
- `BoundsError`: use firm ids in `1:econ.n`, including every entry of
  `firm_indices`.
- `combine`: choose `:min`, `:upstream`, or `:downstream`.
- `components`: choose `:none`, `:upstream`, `:downstream`, or `:both`.
- `shock`: provide `econ.n` finite values in `[0,1]`.
- `final_weights`: provide `econ.n` finite values `>= 0`.
- With `firm_indices`, zero values fill the remaining entries by design.
- In `esri(econ, firm_idx; shock=psi)`, `psi` defines the complete scenario and
  `firm_idx` selects the single-scenario overload. The default closure uses
  `firm_idx` when the call uses the default scenario.

For benchmark scope and the submitted runtime snapshot, see
[Performance](performance.md).
