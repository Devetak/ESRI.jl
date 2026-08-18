# ESRIcascade.jl

[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://devetak.github.io/ESRIcascade.jl/)

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
```

## Essential and non-essential inputs

The one-argument `IndustryInfo(industry_of_firm)` uses a purely linear production
function. To add essential inputs, choose a classification, then use the same
`econ` and `esri` call:

```julia
# 1. Old essential/non-essential behavior, using the new matrix API.
industry_of_firm = rand(1:4, N)
essential_industry = [true, true, false, false]
input_classification = repeat(
    UInt8.(essential_industry) .+ UInt8(1),
    1,
    length(essential_industry),
)
info = IndustryInfo(industry_of_firm, input_classification)

# 2. Bundled IHS matrix: no file loading required.
# codes = ihs_industry_codes() # ["0111", "0112", ..., "9999"]
# industry_id = Dict(code => i for (i, code) in pairs(codes)) # "0111" => 1, "0112" => 2, ...
# industry_of_firm = [industry_id[code] for code in firm_industry_codes] # e.g. ["0111", "0112", ...] -> [1, 2, ...]
# info = IndustryInfo(industry_of_firm, ihs_input_classification()) # use the bundled 616 x 616 matrix

# 3. Custom matrix: rows = suppliers, columns = customers.
# input_classification = [2 1 0; 0 2 1; 1 0 2]
# industry_of_firm = rand(1:3, N)
# info = IndustryInfo(industry_of_firm, input_classification)

econ = ESRIEconomy(W, info)
scores = esri(econ)
```

In the production function:

- `0`: not relevant for production.
- `1`: non-essential input in the linear part of the production function.
- `2`: essential input in the Leontief part of the production function.

Example score distribution from the same kind of run:

![Histogram of example ESRI scores](docs/src/assets/scores_hist.svg)

If most firms are near zero, most single-firm failures have limited economy-wide spillovers. If the histogram has a heavier right tail, some firm failures create much larger losses across the economy.

Build `ESRIEconomy` once and reuse it on the same network.

## Key calls

- `esri(econ; ...)` computes the default single-firm closure for each selected firm. If `firm_indices` is set, unrequested entries stay zero.
- `esri(econ, firm_idx; ...)` solves one scenario and returns a scalar, a named tuple, or `ESRIResult`.
- `esri_shock(econ, psi; ...)` solves one scenario from an explicit capacity cap vector `psi ∈ [0,1]^N`.
- `details=true` or `components=:upstream|:downstream|:both` returns the converged component states for a single scenario.
- `final_weights` replaces the weights in the final ESRI reduction. With custom weights, the scalar is the weighted loss divided by `sum(final_weights)`.
- `shock=psi` on `esri(econ, firm_idx; ...)` replaces the default closure.
- `compute_esri(...)` and `compute_esri_shock(...)` are matrix-first wrappers that build `ESRIEconomy` and dispatch to the matching call.

## Full ESRI C++ reference benchmark

`benchmark/` contains a one-to-one prepared sparse comparison against
`fastcascade::GL_cascade_dynamics_cpp`. It closes every firm at four network
sizes under IHS, legacy essential/non-essential, and all-linear input
classifications. The Julia runner checks every score and the total ESRI against
the C++ result.
See [benchmark/README.md](benchmark/README.md) for the commands and prerequisites.

### Submitted runtime snapshot

The table below is the supplied benchmark snapshot, recorded here without
rerunning the benchmark. Each size uses the same deterministic directed
truncated power-law network specification (mean out-degree 8, exponent 2.3,
degree cap 128, seed 42, no self-links), and reports prepared solve time per
firm. Speedup is
`fastcascade` divided by `ESRIcascade.jl`.
The supplied run uses the benchmark runner's default one-thread BLAS/OpenMP
restriction.

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
tolerance. The detailed benchmark scope and rerun commands are in
[`benchmark/README.md`](benchmark/README.md).

## Reference

Diem, C. et al. *Quantifying firm-level economic systemic risk from nation-wide supply networks* (Scientific Reports, 2022): https://www.nature.com/articles/s41598-022-11522-z

## License

ESRIcascade.jl code is released under the MIT License. The bundled IHS matrix
and labels are separately CC BY 4.0; see `data/LICENSE.md` for attribution.
