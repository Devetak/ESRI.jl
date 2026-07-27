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

Choose one input classification, then use the same `econ` and `esri` call:

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
# codes = ihs_industry_codes()
# industry_id = Dict(code => i for (i, code) in pairs(codes))
# industry_of_firm = [industry_id[code] for code in firm_industry_codes]
# info = IndustryInfo(industry_of_firm, ihs_input_classification())

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
- `final_weights` replaces the weights in the final ESRI reduction. With custom weights, the scalar is the weighted loss divided by `sum(final_weights)`.
- `shock=psi` on `esri(econ, firm_idx; ...)` replaces the default closure. 

## Sparse C++ reference benchmark

`benchmark/` contains a one-to-one prepared sparse comparison against
`fastcascade::GL_cascade_dynamics_cpp`: the same IHS classification matrix,
1,232 firms, 9,856 sparse links, 16 closures, tolerance, and one-thread setting.
Both runners check the same C++ reference scores before printing medians. See
[benchmark/README.md](benchmark/README.md) for the two commands and prerequisites.

## Reference

Diem, C. et al. *Quantifying firm-level economic systemic risk from nation-wide supply networks* (Scientific Reports, 2022): https://www.nature.com/articles/s41598-022-11522-z

## License

ESRIcascade.jl code is released under the MIT License. The bundled IHS matrix
and labels are separately CC BY 4.0; see `data/LICENSE.md` for attribution.
