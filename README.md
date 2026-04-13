# ESRI.jl

`ESRI.jl` computes firm-level economic systemic risk on directed supply networks.

The package follows the ESRI setup of Diem et al. (2022) with a narrower input contract:

- `W` is a square firm-to-firm weight matrix.
- `IndustryInfo` uses one industry id per firm and one Boolean essentiality flag per industry.
- `psi` is always a capacity-cap vector in `[0,1]^N`.

Build `ESRIEconomy(W, info)` once and reuse it when the network is fixed.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/Devetak/ESRI.jl")
```

For local development:

```julia
using Pkg
Pkg.develop(path = "/path/to/ESRI.jl")
```

## Quick start

```julia
using ESRI, SparseArrays
using LinearAlgebra: I

N = 1_000
W = sprand(N, N, 0.01) + 0.1I
info = IndustryInfo(rand(1:5, N), [true, true, true, false, false])

econ = ESRIEconomy(W, info)
scores = esri(econ; maxiter = 50, tol = 1e-3, threads = true)
value = esri(econ, 10; maxiter = 50, tol = 1e-3)
details = esri(econ, 10; maxiter = 50, tol = 1e-3, details = true)
```

## Key calls

- `esri(econ; ...)` computes the default single-firm closure for each selected firm. If `firm_indices` is set, unrequested entries stay zero.
- `esri(econ, firm_idx; ...)` solves one scenario and returns a scalar, a named tuple, or `ESRIResult`.
- `esri_shock(econ, psi; ...)` solves one scenario from an explicit capacity cap vector `psi ∈ [0,1]^N`.
- `final_weights` changes only the numerator of the final ESRI reduction.
- `shock=psi` on `esri(econ, firm_idx; ...)` replaces the default closure. It does not add a second shock on top.

## Docs

- Source docs: `docs/src/`
- Build locally:
  - `julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'`
  - `julia --project=docs docs/make.jl`

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
julia --project benchmark/sparse_powerlaw_esri.jl 10000 truncated_tail
julia --project benchmark/sparse_powerlaw_esri.jl 10000 heavy_tail
```

Quick full-solve smoke benchmark:

```bash
julia --project test/perf_full_powerlaw_esri.jl
```

Larger manual references:

```bash
julia --project benchmark/sparse_powerlaw_esri.jl 50000 truncated_tail
julia --project benchmark/sparse_powerlaw_esri.jl 50000 heavy_tail
```

These are local reference numbers, not CI guarantees.

## Reference

Diem, C. et al. *Quantifying firm-level economic systemic risk from nation-wide supply networks* (Scientific Reports, 2022): https://www.nature.com/articles/s41598-022-11522-z

## License

MIT. See `LICENSE`.
