# Performance Notes

## Sparse vs Dense

`ESRI.jl` supports both dense and sparse matrices. Sparse matrices are preferred for large supply networks.

- Upstream propagation has dedicated sparse and dense kernels.
- Downstream impact matrices are precomputed once in `ESRIEconomy`.
- Memory usage is typically dominated by impact matrices and per-thread work buffers.

## Threading

Use `threads = true` in economy-wide calls:

```julia
scores = esri(econ; threads = true)
```

Threading parallelizes across firm shocks. For deterministic validation, run with `threads = false`.

## Determinism and Reproducibility

- Set `Random.seed!(...)` before synthetic network generation.
- Keep BLAS single-threaded for dense/sparse numerical comparison tests.
- Use fixed `maxiter` and `tol` in benchmarking and regression tests.

## Reuse `ESRIEconomy`

For repeated scenario analysis, build `ESRIEconomy` once and reuse it:

```julia
econ = ESRIEconomy(W, info)
scores = esri(econ)
single = esri(econ, 10)
scenario = esri_shock(econ, psi)
```

This avoids repeatedly rebuilding normalized impact operators.

## Practical Guidance

- For repeated scenario work, build one `ESRIEconomy` and call `esri`/`esri_shock` many times.
- For very sparse national supply graphs, keep input as `SparseMatrixCSC`.
- Use `firm_indices` when only a subset of firms is needed.
