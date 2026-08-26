# ESRI Julia/C++ benchmarks

The full comparison harness closes every firm in the same deterministic directed power-law network at
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

The full harness uses eight workers by default. Network construction, operator setup,
worker startup, and validation are outside the prepared-solve timing. The R
runner checkpoints completed waves; the Julia runner then compares every score
and the total ESRI before writing the comparison diagnostics.

The R script requires `GLcascade` and `fastcascade` 0.9.3.1.
The benchmark runner restricts BLAS/OpenMP to one thread by default; set
`ESRI_UNRESTRICT_BLAS=1` when an unrestricted BLAS comparison is intentional.

```bash
ESRI_IHS_MATRIX=/path/to/EssMatIHS.csv \
  benchmark/run_full_esri_comparison.sh
```

Results are written under `results/full_esri_matrix_comparison/`.

## Submitted runtime snapshot

The following fixed snapshot is the supplied benchmark result. It is documented
verbatim here; the benchmark was not rerun for this documentation update. Each
size uses a deterministic directed truncated power-law network with mean
out-degree 8, exponent 2.3, degree cap 128, seed 42, and no self-links. Values
are prepared solve milliseconds per firm, and speedup is `fastcascade` divided by
`ESRIcascade.jl`. The supplied run uses the runner's default one-thread
BLAS/OpenMP restriction.

| Size | Production function | Workers | ESRIcascade.jl (ms per firm) | fastcascade (ms per firm) | Speedup |
| --- | --- | ---: | ---: | ---: | ---: |
| 10k | IHS | 1 | 0.777 | 172.422 | 221.91× |
| 10k | IHS | 4 | 0.296 | 61.432 | 207.43× |
| 10k | Legacy | 1 | 1.103 | 210.009 | 190.36× |
| 10k | Legacy | 4 | 0.343 | 72.083 | 210.35× |
| 10k | Linear | 1 | 0.826 | 250.429 | 303.30× |
| 10k | Linear | 4 | 0.293 | 86.772 | 295.69× |
| 100k | IHS | 1 | 9.885 | 1900.390 | 192.25× |
| 100k | IHS | 4 | 4.222 | 587.493 | 139.16× |

The table covers all three production functions at 10,000 firms and IHS at
100,000 firms, each with one and four workers. The Julia and `fastcascade`
scores matched within the benchmark convergence tolerance. The full harness
and its prerequisites remain documented above.

For the overnight power-law matrix (10,000 firms under all three production
specifications and 100,000 firms under IHS, each with 1 and 4 workers), run:

```bash
ESRI_IHS_MATRIX=/path/to/EssMatIHS.csv \
  caffeinate -i benchmark/run_overnight_powerlaw.sh
```

The combined result is written under `results/overnight_powerlaw/`. Both
implementations reuse the same generated network file at each firm count.
