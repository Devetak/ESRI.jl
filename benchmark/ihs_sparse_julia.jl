using DelimitedFiles
using ESRIcascade
using SparseArrays
using Statistics

matrix_path = get(ENV, "ESRI_IHS_MATRIX", "")
isfile(matrix_path) || error("Set ESRI_IHS_MATRIX to the labeled EssMatIHS.csv file")

raw = readdlm(matrix_path, Char(44), String)
labels = raw[2:end, 1]
labels == raw[1, 2:end] || error("IHS row and column labels must have the same order")
classification = parse.(Int, raw[2:end, 2:end])

industries = size(classification, 1)
n = parse(Int, get(ENV, "ESRI_BENCHMARK_FIRMS", string(2 * industries)))
n >= industries || error("ESRI_BENCHMARK_FIRMS must be at least the number of industries")
firms_per_industry, extra_firms = divrem(n, industries)
industry_of_firm = vcat(repeat(1:industries, inner = firms_per_industry), 1:extra_firms)
closure_industries = round.(Int, range(1, industries; length = 16))
closure_firms = [(industry - 1) * firms_per_industry + 1 for industry in closure_industries]
offsets = (1, 37, 91, 173, 311, 509, 701, 997)
supplier = repeat(1:n, inner = length(offsets))
customer = [mod1(i + offset, n) for i = 1:n for offset in offsets]
weights = Float64.(1 .+ mod.(13 .* supplier .+ 7 .* customer, 17))
econ = ESRIEconomy(
    sparse(supplier, customer, weights, n, n),
    IndustryInfo(industry_of_firm, classification),
)

function run_julia(econ)
    return esri(
        econ;
        firm_indices = closure_firms,
        maxiter = 300,
        tol = 1e-2,
        threads = false,
        combine = :min,
    )[closure_firms]
end

default_cpp_scores = [
    0.87837245554794330,
    0.83928173885048851,
    0.84026525189212131,
    0.010491619358180508,
    0.83938875684430192,
    0.83919913726570627,
    0.83871134521356561,
    0.83862604469123081,
    0.83842687824673123,
    0.14782516685618610,
    0.013224577635286240,
    0.0064115243633385468,
    0.0046717737458951943,
    0.83930463691132073,
    0.011700530446422091,
    0.0086842594773570356,
]

cpp_scores = if n == 2 * industries
    default_cpp_scores
else
    reference_path = get(ENV, "ESRI_REFERENCE_SCORES", "")
    isfile(reference_path) || error(
        "Set ESRI_REFERENCE_SCORES to the C++ score file for non-default firm counts",
    )
    vec(readdlm(reference_path, ',', Float64))
end
length(cpp_scores) == 16 || error("C++ reference must contain 16 closure scores")

scores = run_julia(econ)
samples = [
    @elapsed(run_julia(econ)) for _ = 1:parse(Int, get(ENV, "ESRI_BENCHMARK_SAMPLES", "5"))
]
println(
    "firms=$n industries=$industries edges=$(length(weights)) closures=16 closure_firms=$(join(closure_firms, ','))",
)
println(
    "downstream_nnz=$(nnz(econ.downstream_impact_essential) + nnz(econ.downstream_impact_nonessential))",
)
println("julia_prepared_median_s=$(median(samples))")
println("julia_allocated_bytes=$(@allocated run_julia(econ))")
println("max_abs_error_to_cpp=$(maximum(abs.(scores .- cpp_scores)))")
maximum(abs.(scores .- cpp_scores)) <= 1e-12 || error("Julia/C++ parity failed")
