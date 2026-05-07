using DelimitedFiles
using LinearAlgebra

BLAS.set_num_threads(1)
ENV["OMP_NUM_THREADS"] = "1"
ENV["OPENBLAS_NUM_THREADS"] = "1"
ENV["MKL_NUM_THREADS"] = "1"
ENV["VECLIB_MAXIMUM_THREADS"] = "1"

if length(ARGS) != 5
    error("usage: julia economic_systemic_risk_runner.jl <legacy_root> <input_csv> <psi_csv> <maxiter> <output_prefix>")
end

legacy_root = ARGS[1]
input_csv = ARGS[2]
psi_csv = ARGS[3]
maxiter = parse(Int, ARGS[4])
output_prefix = ARGS[5]

include(joinpath(legacy_root, "extern.jl"))
include(joinpath(legacy_root, "functions.jl"))

M = initializeMarket(input_csv)
A = buildArrays(M)
psi_mat = parsePsiMat(psi_csv, M)
parsed_args = Dict(
    "input" => input_csv,
    "output" => output_prefix * "_output.csv",
    "tmax" => maxiter,
    "psi_mat" => psi_csv,
    "timeseries" => false,
)

esri_df = nothing
solve_s = @elapsed begin
    esri_df = ESRI(M, A, psi_mat, parsed_args)
end

scores = Float64.(esri_df.esri)
iters = Int.(esri_df.t)
metrics = [
    ("legacy_total_s", solve_s),
    ("legacy_score_min", minimum(scores)),
    ("legacy_score_max", maximum(scores)),
    ("legacy_max_iter", maximum(iters)),
    ("legacy_mean_iter", sum(iters) / length(iters)),
    ("scenario_count", length(scores)),
]

writedlm(output_prefix * "_metrics.tsv", metrics, '\t')
writedlm(output_prefix * "_scores.tsv", scores, '\t')
