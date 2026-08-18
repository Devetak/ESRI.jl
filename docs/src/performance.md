# Performance

Build `ESRIEconomy(W, info)` once and reuse it. That is the main performance rule.

`ESRIcascade.jl` supports both dense and sparse matrices. Sparse matrices are preferred for large supply networks.

`threads=true` applies only to economy-wide `esri(econ; ...)` and `compute_esri(W, info; ...)`. It parallelizes across shocked firms, not within one single-scenario solve.

## Submitted runtime snapshot

The table below is the supplied benchmark snapshot, documented without
rerunning it. Each size uses a deterministic directed truncated power-law
network with mean out-degree 8, exponent 2.3, degree cap 128, seed 42, and no
self-links. Values are prepared solve milliseconds per firm; speedup is
`fastcascade` divided by `ESRIcascade.jl`. The supplied run uses the benchmark
runner's default one-thread BLAS/OpenMP restriction.

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

The Julia and `fastcascade` scores matched within the benchmark convergence
tolerance. See `benchmark/README.md` for the complete harness description and
commands.
