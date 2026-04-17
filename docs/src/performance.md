# Performance

Build `ESRIEconomy(W, info)` once and reuse it. That is the main performance rule.

`ESRIcascade.jl` supports both dense and sparse matrices. Sparse matrices are preferred for large supply networks.

`threads=true` applies only to economy-wide `esri(econ; ...)` and `compute_esri(W, info; ...)`. It parallelizes across shocked firms, not within one single-scenario solve.
