# Performance

The main cost split is fixed. `ESRIEconomy(W, info)` builds the normalized upstream and downstream operators once. Each subsequent call solves one or more fixed-point problems on that precomputed state. If the network does not change, build `ESRIEconomy` once and reuse it.

For large networks, pass `W` as `SparseMatrixCSC`. The package has sparse kernels for the upstream step and sparse downstream operator storage. Dense input is supported, but it is not the intended mode for large national supply graphs.

`threads = true` applies to economy-wide `esri(econ; ...)` and `compute_esri(W, info; ...)`. It parallelizes across shocked firms. It does not parallelize one single-scenario solve such as `esri(econ, firm_idx; ...)` or `esri_shock(econ, psi; ...)`.

The memory cost in threaded economy-wide runs is higher because each worker keeps its own buffers for upstream state, downstream state, industry sums, and shock vector.

`firm_indices` is the main way to reduce work when only part of the economy-wide shock set is needed.

`verbose = true` shows progress in serial economy-wide runs and iteration logs in single-scenario runs. In threaded economy-wide runs the progress display is disabled and the package ignores `verbose = true`.
