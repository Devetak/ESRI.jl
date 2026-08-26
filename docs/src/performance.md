# Performance

Build `ESRIEconomy(W, info)` once and reuse it. That is the main performance rule.

`ESRIcascade.jl` supports dense and sparse matrices. Sparse matrices suit large supply networks.

Start Julia with more than one thread:

```sh
julia --threads=auto --project
```

Then use `threads=true` with economy-wide `esri(econ; ...)` or
`compute_esri(W, info; ...)`. It uses the available Julia threads across
shocked firms. Each single-scenario solve uses one thread.

## Submitted runtime snapshot

The table below records the supplied benchmark snapshot. Each size uses a
deterministic directed truncated power-law network with mean out-degree 8,
exponent 2.3, degree cap 128, seed 42, and edges between distinct firms. Values
are prepared solve milliseconds per firm; speedup is
`fastcascade` divided by `ESRIcascade.jl`. The supplied run uses the benchmark
runner's default one-thread BLAS/OpenMP setting.

The supplied artifacts identify the run as `20260813T134039Z` on Linux
`7.0.0-28-generic`. The commit, processor, Julia version, R version, and R
package versions remain unrecorded. Use the values as a record of this supplied
run. Expect other systems to differ.

| Size | Production function | Julia threads / R workers | ESRIcascade.jl (ms per firm) | fastcascade (ms per firm) | Speedup |
| --- | --- | ---: | ---: | ---: | ---: |
| 10k | IHS | 1 | 0.777 | 172.422 | 221.91× |
| 10k | IHS | 4 | 0.296 | 61.432 | 207.43× |
| 10k | Legacy | 1 | 1.103 | 210.009 | 190.36× |
| 10k | Legacy | 4 | 0.343 | 72.083 | 210.35× |
| 10k | Linear | 1 | 0.826 | 250.429 | 303.30× |
| 10k | Linear | 4 | 0.293 | 86.772 | 295.69× |
| 100k | IHS | 1 | 9.885 | 1900.390 | 192.25× |
| 100k | IHS | 4 | 4.222 | 587.493 | 139.16× |

`Legacy` is the Boolean essential-industry classification.

The Julia and `fastcascade` scores matched within the benchmark convergence
tolerance of `1e-2`. See the
[benchmark guide](https://github.com/Devetak/ESRIcascade.jl/blob/4e44dfb2c319562792dd90fdf4e5acbf6a3acffd/benchmark/README.md)
for the complete harness description and commands.
