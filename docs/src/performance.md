# Performance

Build `ESRIEconomy(W, info)` once and reuse it. That is the main performance rule.

For large networks, pass `W` as `SparseMatrixCSC`. The sparse kernels are the intended path for large supply graphs.

`threads=true` applies only to economy-wide `esri(econ; ...)` and `compute_esri(W, info; ...)`. It parallelizes across shocked firms, not within one single-scenario solve.

The cost of an economy-wide run still grows with the number of shocked firms because the package solves one fixed point per scenario. Use `firm_indices` when only a subset is needed.

Threaded economy-wide runs use more memory because each worker keeps its own buffers. `verbose=true` shows progress in serial economy-wide runs and iteration logs in single-scenario runs. In threaded economy-wide runs it is ignored.
