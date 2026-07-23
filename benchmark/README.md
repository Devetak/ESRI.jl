# Sparse Julia/C++ comparison

These paired scripts time an equivalent prepared sparse workload: 1,232 firms,
the 616 by 616 IHS `0/1/2` classification matrix, 9,856 directed links, and
the first 16 single-firm closures. The rows of `W` are suppliers and its columns
are customers. Both scripts keep network/classification parsing and sparse
operator construction outside the timed region. The C++ setup drops explicit
zero downstream entries before timing, matching Julia's stored sparse operators.
The C++ routine accepts the 16 closures in one native call; Julia's public serial
API completes the same 16 closures with a reused workspace.

Set `ESRI_IHS_MATRIX` to the labeled `EssMatIHS.csv` file, then run each script
with one computational thread:

```bash
ESRI_IHS_MATRIX=/path/to/EssMatIHS.csv \
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
julia --project=. benchmark/ihs_sparse_julia.jl

ESRI_IHS_MATRIX=/path/to/EssMatIHS.csv \
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
Rscript benchmark/ihs_sparse_fastcascade.R
```

The R script requires `GLcascade` and `fastcascade` 0.9.3.1 and calls
`fastcascade::GL_cascade_dynamics_cpp` directly. The Julia script calls only
the public `IndustryInfo`, `ESRIEconomy`, and `esri` API. Both check the same
16 C++ reference scores to `1e-12`; compare their printed medians only after
both parity checks pass.

For the paper-scale firm counts `10_000`, `20_000`, `50_000`, and `100_000`,
set `ESRI_BENCHMARK_FIRMS`. Firms are assigned evenly across the same 616
industries, with any remainder assigned once to the first industries. The C++
script writes the 16 reference scores and Julia reads them before timing:

```bash
for firms in 10000 20000 50000 100000; do
  reference=/tmp/esri-cpp-${firms}.csv
  ESRI_BENCHMARK_FIRMS=$firms ESRI_REFERENCE_SCORES=$reference \
    ESRI_IHS_MATRIX=/path/to/EssMatIHS.csv Rscript benchmark/ihs_sparse_fastcascade.R
  ESRI_BENCHMARK_FIRMS=$firms ESRI_REFERENCE_SCORES=$reference \
    ESRI_IHS_MATRIX=/path/to/EssMatIHS.csv JULIA_NUM_THREADS=1 \
    julia --project=. benchmark/ihs_sparse_julia.jl
done
```
