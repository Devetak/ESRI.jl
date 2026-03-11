# ESRI.jl

Efficient computation of the Economic Systemic Risk Index (ESRI) for large production networks. ESRI.jl is an open-source Julia package that provides a high-level API to compute systemic risk on large, sparse production networks.

At the core of the package is a single high-level API:

- `compute_esri(weight_matrix, info; maxiter=100, tol=1e-2, verbose=false)`  
  where:
  - `weight_matrix` is a dense `Matrix` or sparse `SparseMatrixCSC` with supplier weights, and
  - `info` is an `IndustryInfo` describing industry membership and which industries are considered essential.

`IndustryInfo` encodes, for each firm, an industry identifier and a Boolean flag for whether that industry is essential. `compute_esri` then iteratively propagates shocks through the production network until it converges to an ESRI vector.

## Installation

Until ESRI.jl is registered in the General registry, you can install it directly from the Git repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/Devetak/ESRI.jl")
```

Once registered, you will be able to install it with:

```julia
using Pkg
Pkg.add("ESRI")
```

To work in a local checkout (recommended for development and benchmarking):

```julia
using Pkg
Pkg.develop(path = "/path/to/ESRI.jl")
```

## Quick start

```julia
using ESRI, SparseArrays

N = 1_000
W = sprand(N, N, 0.01) + 0.1I # supplier weights
types = rand(1:5, N)
essential = [true, true, true, false, false]
info = IndustryInfo(types, essential)

esri = compute_esri(W, info; tol=1e-3, maxiter=50)
```

## Public API overview

The main public entry points are:

- `IndustryInfo(industry_ids::AbstractVector{<:Integer}, essential_industry::AbstractVector{Bool})`  
  Construct metadata describing industry membership and which industries are essential.

- `compute_esri(weight_matrix, info; maxiter=100, tol=1e-2, verbose=false)`  
  Compute the Economic Systemic Risk Index for all firms given a (dense or sparse) supplier weight matrix and an `IndustryInfo` instance.

Internally, the implementation is organized in:

- `src/ESRI.jl`: module entry point and exports,
- `src/types.jl`: definitions of `IndustryInfo` and related types,
- `src/esri_computation.jl`: core ESRI iterative computation,
- `src/downstream.jl`, `src/upstream.jl`, `src/impact.jl`: helpers for directional impacts and aggregation.

## Project setup

Instantiate the project environment:

```julia
julia --project -e 'using Pkg; Pkg.instantiate()'
```


## License

ESRI.jl is open source and released under the MIT License. See `LICENSE` for details.

