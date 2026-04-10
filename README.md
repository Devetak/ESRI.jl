# ESRI.jl

ESRI.jl computes, for each firm, the share of the economy that depends on that firm in `[0, 1]`.

## Documentation

- Source docs (with equations): `docs/src/`
- Build locally:
  - `julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'`
  - `julia --project=docs docs/make.jl`

## Installation

Install from this repository (current work-in-progress):

```julia
using Pkg
Pkg.add(url = "https://github.com/Devetak/ESRI.jl")
```

After registration in the General registry:

```julia
using Pkg
Pkg.add("ESRI")
```

For local development:

```julia
using Pkg
Pkg.develop(path = "/path/to/ESRI.jl")
```

## Quick start (sparse only)

```julia
using ESRI, SparseArrays
using LinearAlgebra: I

N = 1_000
W = sprand(N, N, 0.01) + 0.1I # supplier weights (sparse)
industry_ids = rand(1:5, N)
essential_industry = [true, true, true, false, false]
info = IndustryInfo(industry_ids, essential_industry)

econ = ESRIEconomy(W, info)
esri_all = esri(econ; maxiter = 50, tol = 1e-3, threads = true)
esri_i = esri(econ, 10; maxiter = 50, tol = 1e-3)
details_i = esri(econ, 10; maxiter = 50, tol = 1e-3, details = true)
esri_wrapper = compute_esri(W, info; maxiter = 50, tol = 1e-3)
```

## Feature Flags (new)

```julia
# custom final weights (default: econ.row_sums)
scores_w = esri(econ; final_weights = rand(N))

# final reduction mode (default: :min)
scores_u = esri(econ; combine = :upstream)
scores_d = esri(econ; combine = :downstream)

# single-firm + custom shock vector in [0,1]
psi = ones(N); psi[1:10] .= 0.5
one_firm = esri(econ, 7; shock = psi, details = true)

# direct shock scenario API (no firm index)
scenario = esri_shock(econ, psi; combine = :min)
```

## Performance / execution model

ESRI iterates upstream/downstream dynamics per firm, then aggregates a firm ESRI value.
Use `threads=true` in `esri(econ; ...)` (or `compute_esri`) to parallelize across firms.

## Multiprocessing transparency

Threading is explicit in `esri(econ; ...)`; dense linear algebra may still use BLAS threads.

## Reproducibility

```julia
using Random
Random.seed!(42)
```

## Reference paper

Diem, C. et al. *Quantifying firm-level economic systemic risk from nation-wide supply networks* (Scientific Reports, 2022): https://www.nature.com/articles/s41598-022-11522-z

## Public API

- `IndustryInfo(industry_ids::AbstractVector{<:Integer}, essential_industry::AbstractVector{Bool})`
- `ESRIEconomy(weight_matrix, info::IndustryInfo)`
- `esri(econ::ESRIEconomy; maxiter=100, tol=1e-2, verbose=false, threads=false, firm_indices=nothing, final_weights=nothing, combine=:min)`
- `esri(econ::ESRIEconomy, firm_idx::Integer; maxiter=100, tol=1e-2, verbose=false, details=false, components=:none, final_weights=nothing, combine=:min, shock=nothing)`
- `esri_shock(econ::ESRIEconomy, shock; maxiter=100, tol=1e-2, verbose=false, details=false, components=:none, final_weights=nothing, combine=:min)`
- `compute_esri(weight_matrix, info; maxiter=100, tol=1e-2, verbose=false, kwargs...)`
- `compute_esri_shock(weight_matrix, info, shock; maxiter=100, tol=1e-2, verbose=false, kwargs...)`

## License

ESRI.jl is open source and released under the MIT License. See `LICENSE` for details.
