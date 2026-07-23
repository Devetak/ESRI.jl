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
industry_of_firm = vcat(
    repeat(1:industries, inner = firms_per_industry),
    1:extra_firms,
)
offsets = (1, 37, 91, 173, 311, 509, 701, 997)
supplier = repeat(1:n, inner = length(offsets))
customer = [mod1(i + offset, n) for i in 1:n for offset in offsets]
weights = Float64.(1 .+ mod.(13 .* supplier .+ 7 .* customer, 17))
econ = ESRIEconomy(
    sparse(supplier, customer, weights, n, n),
    IndustryInfo(industry_of_firm, classification),
)

function run_julia(econ)
    return esri(
        econ;
        firm_indices = 1:16,
        maxiter = 300,
        tol = 1e-2,
        threads = false,
        combine = :min,
    )[1:16]
end

default_cpp_scores = [
    0.87837245554794330, 0.87825862250992404, 0.87664923655232285, 0.87566688889824862,
    0.87524585194239457, 0.87433651919608424, 0.86060405584903987, 0.84701627242865141,
    0.85599031156433647, 0.85741739269669437, 0.85466742541113849, 0.85393912577670528,
    0.84935069451138689, 0.77626516574715354, 0.84447307983109421, 0.84456652158172307,
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
    @elapsed(run_julia(econ)) for _ in 1:parse(Int, get(ENV, "ESRI_BENCHMARK_SAMPLES", "5"))
]
println("firms=$n industries=$industries edges=$(length(weights)) closures=16")
println("downstream_nnz=$(nnz(econ.downstream_impact_essential) + nnz(econ.downstream_impact_nonessential))")
println("julia_prepared_median_s=$(median(samples))")
println("julia_allocated_bytes=$(@allocated run_julia(econ))")
println("max_abs_error_to_cpp=$(maximum(abs.(scores .- cpp_scores)))")
maximum(abs.(scores .- cpp_scores)) <= 1e-12 || error("Julia/C++ parity failed")
