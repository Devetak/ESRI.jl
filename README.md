# ESRIcascade.jl

[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://devetak.github.io/ESRI.jl/)

ESRIcascade.jl computes, for each firm, the share of the economy that depends on that firm in `[0, 1]`.

`ESRIcascade.jl` is a package for computing the Economic Systemic Risk Index for firms in an economy based on the paper by Diem et al.

In this package, `psi[i]` means how much of its normal capacity firm `i` is allowed to use in the scenario you want to study. `psi[i] = 1.0` means normal operation, `psi[i] = 0.0` means the firm is shut down, and values in between mean the firm can still operate, but only partially. This lets you model shocks such as plant closures, energy shortages, sanctions, or transport disruptions and then measure how those shocks spread through the wider economy.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/Devetak/ESRIcascade.jl")
```

After registration in the General registry:

```julia
using Pkg
Pkg.add("ESRIcascade")
```

For local development:

```julia
using Pkg
Pkg.develop(path = "/path/to/ESRIcascade.jl")
```

## Quick start (sparse only)

```julia
using ESRIcascade, SparseArrays
using LinearAlgebra: I

N = 1_000
W = sprand(N, N, 0.01)
W[1:N+1:end] .= 0
info = IndustryInfo(rand(1:4, N), [true, true, false, false]) # industry 1 and 2 are essential

econ = ESRIEconomy(W, info) # set up the economy
scores = esri(econ; maxiter = 40, tol = 1e-3) # compute ESRI for each firm
nothing
```

Example score distribution from the same kind of run:

![Histogram of example ESRI scores](docs/src/assets/scores_hist.svg)

If most firms are near zero, most single-firm failures have limited economy-wide spillovers. If the histogram has a heavier right tail, some firm failures create much larger losses across the economy.

Build `ESRIEconomy` once and reuse it on the same network.

## Key calls

- `esri(econ; ...)` computes the default single-firm closure for each selected firm. If `firm_indices` is set, unrequested entries stay zero.
- `esri(econ, firm_idx; ...)` solves one scenario and returns a scalar, a named tuple, or `ESRIResult`.
- `esri_shock(econ, psi; ...)` solves one scenario from an explicit capacity cap vector `psi ∈ [0,1]^N`.
- `final_weights` changes only the numerator of the final ESRI reduction.
- `shock=psi` on `esri(econ, firm_idx; ...)` replaces the default closure. It does not add a second shock on top.

## Reference benchmarks

Local reference runs from `2026-04-12` on `Apple M2`, with `JULIA_NUM_THREADS=1`, `mean_degree=7`, `alpha=2.3`, `nindustries=50`, and `maxiter=30`. These timings call full `esri(econ; ...)` over all firms.

| mode | firms | nnz | max_degree | p99_degree | top1pct_edge_share | build_s | solve_s |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `truncated_tail` | 5_000 | 35_000 | 111 | 25 | 0.0637 | 0.0312 | 6.5329 |
| `heavy_tail` | 5_000 | 35_000 | 385 | 26 | 0.0928 | 0.0020 | 5.7411 |
| `truncated_tail` | 10_000 | 70_000 | 128 | 24 | 0.0616 | 0.0388 | 26.6213 |
| `heavy_tail` | 10_000 | 70_000 | 964 | 25 | 0.1071 | 0.0334 | 26.4600 |

Run with:

```bash
julia --project test/perf_full_powerlaw_esri.jl 10000 truncated_tail
julia --project test/perf_full_powerlaw_esri.jl 10000 heavy_tail
```

## Reference

Diem, C. et al. *Quantifying firm-level economic systemic risk from nation-wide supply networks* (Scientific Reports, 2022): https://www.nature.com/articles/s41598-022-11522-z

## License

ESRIcascade.jl is open source and released under the MIT License. See `LICENSE` for details.
