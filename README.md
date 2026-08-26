# ESRIcascade.jl

[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://devetak.github.io/ESRIcascade.jl/)

ESRIcascade.jl computes one Economic Systemic Risk Index (ESRI) score per
firm. Each score measures the economy-wide output loss from a shock to that
firm and lies in `[0, 1]`.

The package implements the ESRI framework from Diem et al. for firm-to-firm supply networks.

`psi[i]` sets the capacity fraction available to firm `i` in a scenario.
`1.0` represents normal capacity. `0.0` represents closure. Values between
these bounds represent partial capacity. Use `psi` to model plant closures,
energy shortages, sanctions, and transport disruptions. ESRI measures the
resulting economy-wide propagation.

## Installation

```julia
using Pkg
Pkg.add("ESRIcascade")
```

Requires Julia 1.12 or later.

For the development version:

```julia
using Pkg
Pkg.add(url = "https://github.com/Devetak/ESRIcascade.jl")
```

For local development:

```julia
using Pkg
Pkg.develop(path = "/path/to/ESRIcascade.jl")
```

## Quick start with a sparse network

```julia
using ESRIcascade, SparseArrays, Random
Random.seed!(42)

N = 1_000
W = sprand(N, N, 0.01)
W[1:N+1:end] .= 0
info = IndustryInfo(rand(1:4, N), [true, true, false, false]) # industry 1 and 2 are essential

econ = ESRIEconomy(W, info) # set up the economy
scores = esri(econ; maxiter = 40, tol = 1e-3) # compute ESRI for each firm
```

`W[i, j]` is the nonnegative supply from supplier `i` to customer `j`.

## Input classifications

The one-argument `IndustryInfo(industry_of_firm)` uses a purely linear production
function. Choose a Boolean vector or square integer-code matrix to classify
inputs. Reuse the same `econ` and `esri` calls:

```julia
# 1. Legacy Boolean mapping expressed through the matrix API.
industry_of_firm = rand(1:4, N)
essential_industry = [true, true, false, false]
input_classification = repeat(
    UInt8.(essential_industry) .+ UInt8(1),
    1,
    length(essential_industry),
)
info = IndustryInfo(industry_of_firm, input_classification)

# 2. Bundled IHS matrix.
# codes = ihs_industry_codes() # ["0111", "0112", ..., "9999"]
# industry_id = Dict(code => i for (i, code) in pairs(codes)) # "0111" => 1, "0112" => 2, ...
# industry_of_firm = [industry_id[code] for code in firm_industry_codes] # e.g. ["0111", "0112", ...] -> [1, 2, ...]
# info = IndustryInfo(industry_of_firm, ihs_input_classification()) # use the bundled 616 × 616 matrix

# 3. Custom matrix: rows = suppliers, columns = customers.
# input_classification = UInt8[2 1 0; 0 2 1; 1 0 2]
# industry_of_firm = rand(1:3, N)
# info = IndustryInfo(industry_of_firm, input_classification)

econ = ESRIEconomy(W, info)
scores = esri(econ)
```

The production function uses these input classes:

- `0`: no direct downstream term; the input remains in `W`, upstream
  propagation, and the class-`1` normalization denominator.
- `1`: linear downstream input.
- `2`: essential downstream input.

Example score distribution from the quick start:

![Histogram of example ESRI scores](docs/src/assets/scores_hist.svg)

Most firms near zero indicate small economy-wide spillovers from most
single-firm shocks. A heavier right tail indicates larger losses from selected
firm shocks.

Build `ESRIEconomy` once and reuse it on the same network.

## Key calls

- `esri(econ; ...)` computes the default single-firm closure for each selected
  firm. Set `firm_indices` to select firms; the remaining entries stay zero.
- `esri(econ, firm_idx; ...)` solves one scenario and returns a scalar, a named
  tuple, or `ESRIResult`.
- `esri_shock(econ, psi; ...)` solves one scenario from an explicit capacity
  cap vector `psi ∈ [0,1]^N`.
- `details=true` or `components=:upstream|:downstream|:both` returns the
  final component states for a single scenario.
- `final_weights` replaces the weights in the final ESRI reduction. With
  custom weights, the scalar is the weighted loss divided by
  `sum(final_weights)`. A zero total weight returns `0`.
- `shock=psi` on `esri(econ, firm_idx; ...)` supplies the complete capacity-cap scenario.
- `compute_esri(...)` and `compute_esri_shock(...)` are matrix-first wrappers.
  They build `ESRIEconomy` and dispatch to the matching call.

## Full ESRI C++ reference benchmark

`benchmark/` contains a one-to-one prepared sparse comparison against
`fastcascade::GL_cascade_dynamics_cpp`. It closes every firm at four network
sizes under IHS, legacy Boolean essentiality, and all-linear input
classifications. The Julia runner checks every score and the total ESRI against
the C++ result.
See [benchmark/README.md](benchmark/README.md) for the commands and prerequisites.

### Submitted runtime snapshot

The table below records the supplied benchmark snapshot. Each size uses the
same deterministic directed truncated power-law network specification. The
network has mean out-degree 8, exponent 2.3,
degree cap 128, seed 42, and edges between distinct firms. It reports prepared solve time per
firm. Speedup is
`fastcascade` divided by `ESRIcascade.jl`.
The supplied run uses the benchmark runner's default one-thread BLAS/OpenMP
setting.

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

Diem, C. et al. [*Quantifying firm-level economic systemic risk from
nation-wide supply networks*](https://www.nature.com/articles/s41598-022-11522-z)
(Scientific Reports, 2022).

## License

ESRIcascade.jl code is released under the MIT License. The bundled IHS matrix
and labels are separately CC BY 4.0; see [data/LICENSE.md](data/LICENSE.md) for
attribution.
