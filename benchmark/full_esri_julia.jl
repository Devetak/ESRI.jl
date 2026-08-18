using DelimitedFiles
using ESRIcascade
using LinearAlgebra
using Printf
using SparseArrays

get(ENV, "ESRI_UNRESTRICT_BLAS", "0") == "1" || BLAS.set_num_threads(1)

classification_name = get(ENV, "ESRI_CLASSIFICATION", "ihs")
classification_name in ("ihs", "legacy", "linear") ||
    error("ESRI_CLASSIFICATION must be ihs, legacy, or linear")
n = parse(Int, get(ENV, "ESRI_BENCHMARK_FIRMS", "10000"))
workers = parse(Int, get(ENV, "ESRI_BENCHMARK_WORKERS", "8"))
convergence_tol = parse(Float64, get(ENV, "ESRI_CONVERGENCE_TOL", "1e-2"))
isfinite(convergence_tol) && convergence_tol > 0 || error("ESRI_CONVERGENCE_TOL must be positive")
workers <= Threads.nthreads() ||
    error("Start Julia with at least $workers threads")
output_root = get(ENV, "ESRI_OUTPUT_DIR", joinpath("results", "full_esri_matrix_comparison"))
network_path = get(ENV, "ESRI_NETWORK_FILE", "")
isfile(network_path) || error("Set ESRI_NETWORK_FILE to the shared power-law edge CSV")
case_dir = joinpath(output_root, "$(classification_name)_$(n)")
mkpath(case_dir)

ihs = ihs_input_classification()
industries = size(ihs, 1)
n >= industries || error("ESRI_BENCHMARK_FIRMS must be at least $industries")

classification = if classification_name == "ihs"
    ihs
elseif classification_name == "legacy"
    essential = iseven.(1:industries)
    repeat(reshape(UInt8.(essential) .+ UInt8(1), :, 1), 1, industries)
else
    fill(UInt8(1), industries, industries)
end

expected_profiles = classification_name == "ihs" ? 56 : 1
length(unique(collect(eachcol(classification)))) == expected_profiles ||
    error("Unexpected number of input-classification profiles")

firms_per_industry, extra_firms = divrem(n, industries)
industry_of_firm = vcat(repeat(1:industries, inner = firms_per_industry), 1:extra_firms)
network_data = readdlm(network_path, ',', Float64)
size(network_data, 2) == 3 || error("Power-law edge CSV must have three columns")
supplier = Int.(network_data[:, 1])
customer = Int.(network_data[:, 2])
weights = network_data[:, 3]
length(supplier) == 8n || error("Expected $(8n) power-law links")
all(x -> 1 <= x <= n, supplier) && all(x -> 1 <= x <= n, customer) ||
    error("Network ID out of range")
all(supplier .!= customer) || error("Power-law network contains self-links")
all(isfinite, weights) && all(weights .> 0) || error("Network weights must be positive")
weight_matrix = sparse(supplier, customer, weights, n, n)
nnz(weight_matrix) == length(weights) || error("Power-law network contains duplicate links")
econ = ESRIEconomy(
    weight_matrix,
    IndustryInfo(industry_of_firm, classification),
)
network_data = supplier = customer = weights = weight_matrix = nothing
GC.gc()

run_part(ids) = esri(
    econ;
    firm_indices = ids,
    maxiter = 300,
    tol = convergence_tol,
    threads = false,
    combine = :min,
)

run_part([1]) # compile outside the timed region
parts = [collect(worker:workers:n) for worker = 1:workers]
GC.gc()
started = time_ns()
tasks = [Threads.@spawn run_part(ids) for ids in parts]
partial_scores = fetch.(tasks)
julia_total_s = (time_ns() - started) / 1e9

scores = zeros(n)
for (ids, values) in zip(parts, partial_scores)
    scores[ids] = values[ids]
end

cpp_path = joinpath(case_dir, "cpp_scores.csv")
isfile(cpp_path) || error("Run full_esri_fastcascade.R first")
cpp = readdlm(cpp_path, ',', Float64)
size(cpp) == (n, 2) || error("C++ score file must contain exactly one row per firm")
cpp_ids = Int.(cpp[:, 1])
cpp_scores = cpp[:, 2]
cpp_ids == collect(1:n) || error("C++ scores are not in firm order")

open(joinpath(case_dir, "julia_scores.csv"), "w") do io
    for i in eachindex(scores)
        @printf(io, "%d,%.17g\n", i, scores[i])
    end
end

score_differences = abs.(scores .- cpp_scores)
max_abs_score_difference = maximum(score_differences)
max_difference_firm = argmax(score_differences)
mean_abs_score_difference = sum(score_differences) / n
within_tolerance = max_abs_score_difference <= convergence_tol
julia_total_esri = sum(scores)
cpp_total_esri = sum(cpp_scores)
abs_total_esri_difference = abs(julia_total_esri - cpp_total_esri)
within_tolerance || @warn "Julia/FastCascade score difference exceeds convergence tolerance" max_abs_score_difference convergence_tol

summary = Dict(
    split(line, '='; limit = 2) for line in readlines(joinpath(case_dir, "cpp_summary.txt"))
)
cpp_total_s = parse(Float64, summary["cpp_total_s"])

open(joinpath(case_dir, "comparison.csv"), "w") do io
    println(
        io,
        "classification,firms,workers,julia_total_s,cpp_total_s,speedup," *
        "julia_total_esri,cpp_total_esri,max_abs_score_difference," *
        "max_difference_firm,julia_score_at_max_difference,cpp_score_at_max_difference," *
        "mean_abs_score_difference,abs_total_esri_difference,comparison_tolerance," *
        "within_tolerance",
    )
    @printf(
        io,
        "%s,%d,%d,%.9f,%.9f,%.9f,%.17g,%.17g,%.17g,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%s\n",
        classification_name,
        n,
        workers,
        julia_total_s,
        cpp_total_s,
        cpp_total_s / julia_total_s,
        julia_total_esri,
        cpp_total_esri,
        max_abs_score_difference,
        max_difference_firm,
        scores[max_difference_firm],
        cpp_scores[max_difference_firm],
        mean_abs_score_difference,
        abs_total_esri_difference,
        convergence_tol,
        string(within_tolerance),
    )
end

println(
    "classification=$classification_name firms=$n workers=$workers " *
    "julia_total_s=$julia_total_s total_esri=$julia_total_esri",
)
