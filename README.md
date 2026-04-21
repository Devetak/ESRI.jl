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

Use the native comparison harness to benchmark `ESRIcascade` and the subset of Diem's bundled `fastcascade` C++ implementation that matches this package's supported setup: one global essential/non-essential sector split, no substitution, and default single-firm shock scenarios on the same sparse power-law network.

| firms | nnz | max_degree | p99_degree | top1pct_edge_share | esri_build_s | esri_solve_s | diem_total_s | esri_total_speedup_x |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `10_000` | run `benchmark/compare_fastcascade.jl` | | | | | | | |
| `50_000` | run `benchmark/compare_fastcascade.jl` | | | | | | | |
| `100_000` | run `benchmark/compare_fastcascade.jl` | | | | | | | |

Run with:

```bash
julia --project=. test/perf_full_powerlaw_esri.jl 1000
julia --project=. benchmark/compare_fastcascade.jl 1000
julia --project=. benchmark/compare_fastcascade.jl 10000 50000 100000
```

The `test/perf_full_powerlaw_esri.jl` command is a quick ESRIcascade-only smoke check. The `benchmark/compare_fastcascade.jl` command generates one shared sparse power-law network per requested size, runs both implementations on that exact input, and prints a Markdown table you can paste back into the README.

## Reference

Diem, C. et al. *Quantifying firm-level economic systemic risk from nation-wide supply networks* (Scientific Reports, 2022): https://www.nature.com/articles/s41598-022-11522-z

## License

ESRIcascade.jl is open source and released under the MIT License. See `LICENSE` for details.
