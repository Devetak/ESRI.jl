# Bundled IHS classification

`ihs_classification.bin` is a compact, lossless representation of a labeled
616 by 616 IHS input-classification matrix. It stores 56 unique
customer-industry profiles plus their column map; `ihs_industry_codes.txt`
stores the row/column order.

The direct provenance and license source is the CC BY 4.0 version 3 Zenodo
replication archive for Pichler et al. (2022), DOI
[`10.5281/zenodo.7113762`](https://doi.org/10.5281/zenodo.7113762). Its
`0_basic_data_IHS_matrices.R` processes the IHS industry-analyst survey and
the archive contains its processed 55 by 55 `A.essential3.csv` output. The
archive checksum published by Zenodo is
`ea69b9a9ae727417353eb9d52397ad39` (MD5); the source entry's SHA-256 is
`fe30f4ccc68493fe8e773846ca2b2134e0e4d11d7355b88c9990bd936c4aab22`.

`verify_zenodo_ihs.jl` expands that file by NACE Rev. 2 division, maps its
`0/0.5/1` values to `0/1/2`, and verifies every 616 by 616 value plus the
bundled label order. Run it after downloading the pinned archive:

```sh
julia --project=. data/verify_zenodo_ihs.jl /path/to/covid19inputoutput.zip
```

The check uses SHA-256
`5f90f03db4af0c144e989bb46d21db016ef7841bf5d9ff4e1c23de79a97f3a27`
for the reconstructed matrix and
`29cbacaf97a3cb51c9e0b92603a024f3b1f9f6a45392c0c675a4bc07142ac61b`
for the newline-delimited labels.

The two bundled data assets are CC BY 4.0; see [LICENSE.md](LICENSE.md) for
the required attribution and modification notice. The package code remains
MIT-licensed.

The package exposes these assets through `ihs_input_classification()` and
`ihs_industry_codes()`.

The bundled IHS classification is opt-in: `IndustryInfo(industry_of_firm)`
uses the purely linear baseline, while passing `ihs_input_classification()`
enables the bundled essential/non-essential classification.

The separate [ICIO regression fixture](ICIO/README.md) is used only by
`test/test_ESRI_ICIO.jl`; it is not loaded by the package at runtime.

For the submitted runtime comparison using this bundled data, see the
[performance benchmark snapshot](../docs/src/performance.md).
