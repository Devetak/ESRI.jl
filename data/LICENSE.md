# License for the bundled IHS classification data

`ihs_classification.bin` and `ihs_industry_codes.txt` are licensed under
[Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/)
(CC BY 4.0). The package code remains licensed under the repository's MIT
License.

## Attribution

The files are derived from the version 3 replication archive for:

> Pichler, Anton; Pangallo, Marco; del Rio-Chanona, R. Maria; Lafond,
> François; and Farmer, J. Doyne (2022). *Forecasting the propagation of
> pandemic shocks with a dynamic input-output model*. Zenodo.
> https://doi.org/10.5281/zenodo.7113762

The Zenodo record identifies the archive as CC BY 4.0. Its
`0_basic_data_IHS_matrices.R` processes the IHS industry-analyst survey; this
package expands its processed `A.essential3.csv` to NACE Rev. 2 four-digit
industries, converts `0/0.5/1` to `0/1/2`, and packs the result for loading.
IHS Markit is identified in the source archive as the origin of the underlying
survey. No endorsement by any source or contributor is implied.
