using DelimitedFiles
using ESRIcascade
using LinearAlgebra
using SparseArrays

BLAS.set_num_threads(1)
ENV["OMP_NUM_THREADS"] = "1"
ENV["OPENBLAS_NUM_THREADS"] = "1"
ENV["MKL_NUM_THREADS"] = "1"
ENV["VECLIB_MAXIMUM_THREADS"] = "1"

function _read_vector(path, ::Type{T}) where {T}
    data = readdlm(path, '\t', T)
    return vec(data)
end

function _read_edges(path)
    rows = readdlm(path, '\t', Any)
    I = Int.(rows[:, 1])
    J = Int.(rows[:, 2])
    V = Float64.(rows[:, 3])
    return I, J, V
end

function main(args = ARGS)
    if length(args) != 6
        error("usage: julia esri_speed_runner.jl <prepared_dir> <tol> <maxiter> <threaded> <metrics_out> <scores_out>")
    end

    prepared_dir = args[1]
    tol = parse(Float64, args[2])
    maxiter = parse(Int, args[3])
    threaded_raw = lowercase(args[4])
    metrics_out = args[5]
    scores_out = args[6]

    threaded =
        threaded_raw in ("true", "1", "yes", "on") ? true :
        threaded_raw in ("false", "0", "no", "off") ? false :
        throw(ArgumentError("threaded must be true or false"))

    industries = Int.(_read_vector(joinpath(prepared_dir, "industries.tsv"), Int))
    row_weights = Float64.(_read_vector(joinpath(prepared_dir, "row_weights.tsv"), Float64))
    essential_flags = Bool.(Int.(_read_vector(joinpath(prepared_dir, "essential.tsv"), Int)))
    I, J, V = _read_edges(joinpath(prepared_dir, "edges.tsv"))

    n = length(industries)
    W = sparse(I, J, V, n, n)
    info = IndustryInfo(industries, essential_flags)

    build_s = @elapsed econ = ESRIEconomy(W, info)
    scores = nothing
    solve_s = @elapsed begin
        scores = esri(
            econ;
            final_weights = row_weights,
            maxiter = maxiter,
            tol = tol,
            threads = threaded,
            combine = :min,
        )
    end

    metrics = [
        ("esri_build_s", build_s),
        ("esri_solve_s", solve_s),
        ("esri_score_min", minimum(scores)),
        ("esri_score_max", maximum(scores)),
        ("esri_nthreads", Threads.nthreads()),
        ("esri_threaded", threaded ? 1.0 : 0.0),
        ("scenario_count", length(scores)),
    ]

    writedlm(metrics_out, metrics, '\t')
    writedlm(scores_out, scores, '\t')
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
