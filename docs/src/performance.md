# Performance

Build `ESRIEconomy(W, info)` once and reuse it. That is the main performance rule.

Pass `W` as `SparseMatrixCSC`. The sparse kernels are the intended path for large supply graphs.

`threads=true` applies only to economy-wide `esri(econ; ...)` and `compute_esri(W, info; ...)`. It parallelizes across shocked firms, not within one single-scenario solve.