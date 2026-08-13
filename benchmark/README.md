# Full ESRI Julia/C++ benchmark

The benchmark closes every firm in the same deterministic directed power-law network at
10,000, 20,000, 50,000, and 100,000 firms. It compares three 616-industry
classifications: the IHS `0/1/2` matrix, the legacy essential/non-essential
rule with alternating essential supplier industries, and an all-linear matrix.

Out-degrees follow a truncated power law with exponent 2.3, mean 8, and cap
128. The generator uses seed 42, excludes self-links, and writes one shared
edge list that both Julia and FastCascade/R load.

Both solvers use the tutorial convergence tolerance of `1e-2`. Because their
stopping criteria differ, the comparison reports maximum and mean absolute
score differences, the firm and two scores at the maximum difference, and a
`within_tolerance` flag instead of requiring equality. Firm IDs must match
one-to-one before the comparison runs.

Both implementations use eight workers. Network construction, operator setup,
worker startup, and validation are outside the prepared-solve timing. The R
runner checkpoints completed waves; the Julia runner then compares every score
and the total ESRI before writing the comparison diagnostics.

The R script requires `GLcascade` and `fastcascade` 0.9.3.1.

```bash
ESRI_IHS_MATRIX=/path/to/EssMatIHS.csv \
  benchmark/run_full_esri_comparison.sh
```

Results are written under `results/full_esri_matrix_comparison/`.

For the overnight power-law matrix (10,000 firms under all three production
specifications and 100,000 firms under IHS, each with 1 and 4 workers), run:

```bash
ESRI_IHS_MATRIX=/path/to/EssMatIHS.csv \
  caffeinate -i benchmark/run_overnight_powerlaw.sh
```

The combined result is written under `results/overnight_powerlaw/`. Both
implementations reuse the same generated network file at each firm count.
