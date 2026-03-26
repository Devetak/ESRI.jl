# ESRI.jl

ESRI.jl computes, for each firm, the share of the economy that depends on that firm in `[0, 1]`.

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
using ESRI
using SparseArrays
using LinearAlgebra: I

N = 1_000
W = sprand(N, N, 0.01) + 0.1I # supplier weights (sparse)

industry_ids = rand(1:5, N)
essential_industry = [true, true, true, false, false] # length = number of industries
info = IndustryInfo(industry_ids, essential_industry)

esri = compute_esri(W, info; maxiter = 50, tol = 1e-3)
```

## Performance / execution model

ESRI is computed by iterating upstream and downstream shock dynamics for each firm, then aggregating an ESRI contribution for that firm. By default it runs serially over firms; you can opt into Julia threading via an API keyword on `compute_esri`.

## Multiprocessing transparency

`compute_esri` controls parallelism explicitly. If you enable threading, firms are computed independently and written to unique output slots; dense linear algebra inside the method may still use BLAS threads depending on your Julia/BLAS setup.

## Reproducibility

```julia
using Random, ESRI
using SparseArrays
using LinearAlgebra: I

Random.seed!(42)

N = 1_000
W = sprand(N, N, 0.01) + 0.1I

industry_ids = rand(1:5, N)
essential_industry = [true, true, true, false, false]
info = IndustryInfo(industry_ids, essential_industry)

esri = compute_esri(W, info; maxiter = 50, tol = 1e-3)
```

## Reference paper

Diem, C. et al. *Quantifying firm-level economic systemic risk from nation-wide supply networks* (Scientific Reports, 2022): https://www.nature.com/articles/s41598-022-11522-z

## Public API

- `IndustryInfo(industry_ids::AbstractVector{<:Integer}, essential_industry::AbstractVector{Bool})`
- `compute_esri(weight_matrix, info; maxiter = 100, tol = 1e-2, verbose = false, kwargs...)`

## License

ESRI.jl is open source and released under the MIT License. See `LICENSE` for details.

