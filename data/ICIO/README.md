# ICIO regression fixture

These committed CSV files provide the reference data for
`test/test_ESRI_ICIO.jl`. They are test fixtures, not part of the package's
runtime data-loading API.

- `p_icio.csv` contains 3,465 production-industry codes.
- `W_ICIO.csv` contains the corresponding 3,465 by 3,465 firm-to-firm weight
  matrix.
- `ESRI_icio.csv` contains three reference score columns, checked in the order
  `:min`, `:downstream`, and `:upstream`.

The test maps the first two digits of each production code to the bundled
`ihs_industry_codes()`, constructs
`IndustryInfo(industry_ids, ihs_input_classification())`, and compares all
three score vectors with an absolute tolerance of `1e-6`.

Thanks to Jan Filankowski for contributing this ICIO reference fixture.
