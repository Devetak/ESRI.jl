using DelimitedFiles
using ESRIcascade
using LinearAlgebra
using SparseArrays

BLAS.set_num_threads(1)
ENV["OMP_NUM_THREADS"] = "1"
ENV["OPENBLAS_NUM_THREADS"] = "1"
ENV["MKL_NUM_THREADS"] = "1"
ENV["VECLIB_MAXIMUM_THREADS"] = "1"

struct TutorialBenchmarkConfig
    scenarios::Int
    tol::Float64
    diff_tol::Float64
    maxiter::Int
    modes::Vector{Symbol}
    gl_ncores::Union{Bool, Int}
    gl_load_balance::Bool
    julia_threaded::Bool
end

function _parse_bool_arg(raw::AbstractString, name::AbstractString)
    value = lowercase(raw)
    if value in ("true", "1", "yes", "on")
        return true
    elseif value in ("false", "0", "no", "off")
        return false
    end
    throw(ArgumentError("$name must be true or false"))
end

function _parse_gl_ncores_arg(raw::AbstractString)
    value = lowercase(raw)
    if value in ("false", "0")
        return false
    end
    ncores = parse(Int, raw)
    ncores > 0 || throw(ArgumentError("gl-ncores must be a positive integer or false"))
    return ncores
end

function parse_args(args)
    scenarios = 500
    tol = 1e-2
    diff_tol = 1e-2
    maxiter = 100
    modes = Symbol[:linear, :leontief]
    gl_ncores = false
    gl_load_balance = false
    julia_threaded = false

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--scenarios"
            i += 1
            scenarios = parse(Int, args[i])
        elseif arg == "--tol"
            i += 1
            tol = parse(Float64, args[i])
        elseif arg == "--diff-tol"
            i += 1
            diff_tol = parse(Float64, args[i])
        elseif arg == "--maxiter"
            i += 1
            maxiter = parse(Int, args[i])
        elseif arg == "--mode"
            i += 1
            modes = Symbol[Symbol(args[i])]
        elseif arg == "--gl-ncores"
            i += 1
            gl_ncores = _parse_gl_ncores_arg(args[i])
        elseif arg == "--gl-load-balance"
            i += 1
            gl_load_balance = _parse_bool_arg(args[i], "gl-load-balance")
        elseif arg == "--julia-threaded"
            i += 1
            julia_threaded = _parse_bool_arg(args[i], "julia-threaded")
        elseif startswith(arg, "--")
            error("unknown argument: $arg")
        else
            scenarios = parse(Int, arg)
        end
        i += 1
    end

    scenarios > 0 || throw(ArgumentError("scenarios must be positive"))
    tol > 0 || throw(ArgumentError("tol must be positive"))
    diff_tol > 0 || throw(ArgumentError("diff-tol must be positive"))
    maxiter > 0 || throw(ArgumentError("maxiter must be positive"))
    all(in((:linear, :leontief)), modes) || throw(ArgumentError("mode must be linear or leontief"))
    julia_threaded || Threads.nthreads() == 1 || throw(
        ArgumentError("launch Julia with --threads 1 or pass --julia-threaded true for multicore ESRIcascade timing")
    )

    return TutorialBenchmarkConfig(
        scenarios,
        tol,
        diff_tol,
        maxiter,
        unique(modes),
        gl_ncores,
        gl_load_balance,
        julia_threaded,
    )
end

function _script_path(name)
    return joinpath(@__DIR__, name)
end

function _tutorial_root()
    return joinpath(dirname(@__DIR__), "esri_tutorial")
end

function _read_vector(path, ::Type{T}) where {T}
    data = readdlm(path, '\t', T)
    return vec(data)
end

function _read_string_vector(path)
    return vec(string.(readdlm(path, '\t', Any)))
end

function _read_metadata(path)
    metadata = Dict{String, String}()
    for row in eachrow(readdlm(path, '\t', Any))
        metadata[string(row[1])] = string(row[2])
    end
    return metadata
end

function _read_edges(path)
    rows = readdlm(path, '\t', Any)
    I = Int.(rows[:, 1])
    J = Int.(rows[:, 2])
    V = Float64.(rows[:, 3])
    return I, J, V
end

function prepare_inputs(tutorial_root)
    isdir(tutorial_root) || error("expected local tutorial checkout at $tutorial_root")

    prepared_dir = mktempdir()
    prepare_script = _script_path("esri_tutorial_prepare.R")
    run(`Rscript $prepare_script $tutorial_root $prepared_dir`)

    metadata = _read_metadata(joinpath(prepared_dir, "metadata.tsv"))
    industries = Int.(_read_vector(joinpath(prepared_dir, "industries.tsv"), Int))
    row_weights = Float64.(_read_vector(joinpath(prepared_dir, "row_weights.tsv"), Float64))
    industry_labels = _read_string_vector(joinpath(prepared_dir, "industry_labels.tsv"))
    I, J, V = _read_edges(joinpath(prepared_dir, "edges.tsv"))
    n = length(industries)
    W = sparse(I, J, V, n, n)

    return (
        prepared_dir = prepared_dir,
        metadata = metadata,
        W = W,
        industries = industries,
        row_weights = row_weights,
        industry_labels = industry_labels,
    )
end

function run_glcascade(
    prepared_dir,
    mode::Symbol,
    tol::Float64,
    scenarios::Int;
    gl_ncores::Union{Bool, Int} = false,
    gl_load_balance::Bool = false,
)
    metrics_path = joinpath(prepared_dir, "gl_$(mode)_metrics.tsv")
    scores_path = joinpath(prepared_dir, "gl_$(mode)_scores.tsv")
    runner_script = _script_path("esri_tutorial_runner.R")
    gl_ncores_arg = gl_ncores === false ? "false" : string(gl_ncores)
    gl_load_balance_arg = gl_load_balance ? "true" : "false"
    run(
        `Rscript $runner_script $prepared_dir $(String(mode)) $(string(tol)) $(string(scenarios)) $metrics_path $scores_path $gl_ncores_arg $gl_load_balance_arg`
    )
    metrics = Dict{String, Float64}()
    for row in eachrow(readdlm(metrics_path, '\t', Any))
        metrics[string(row[1])] = row[2] isa Number ? Float64(row[2]) : parse(Float64, string(row[2]))
    end
    scores = Float64.(_read_vector(scores_path, Float64))
    return metrics, scores
end

function run_julia_mode(
    W,
    industries,
    row_weights,
    mode::Symbol,
    scenarios::Int,
    tol::Float64,
    maxiter::Int;
    threaded::Bool = false,
)
    nindustries = maximum(industries)
    essential = fill(mode == :leontief, nindustries)
    info = IndustryInfo(industries, essential)
    econ = ESRIEconomy(W, info)
    firm_indices = collect(1:min(scenarios, size(W, 1)))
    score_values = nothing
    elapsed = @elapsed begin
        score_values = esri(
            econ;
            firm_indices = firm_indices,
            final_weights = row_weights,
            maxiter = maxiter,
            tol = tol,
            threads = threaded,
            combine = :min,
        )
    end
    return elapsed, score_values[firm_indices]
end

function compare_mode(inputs, config::TutorialBenchmarkConfig, mode::Symbol)
    julia_solve_s, julia_scores = run_julia_mode(
        inputs.W,
        inputs.industries,
        inputs.row_weights,
        mode,
        config.scenarios,
        config.tol,
        config.maxiter;
        threaded = config.julia_threaded,
    )
    gl_metrics, gl_scores = run_glcascade(
        inputs.prepared_dir,
        mode,
        config.tol,
        config.scenarios;
        gl_ncores = config.gl_ncores,
        gl_load_balance = config.gl_load_balance,
    )

    length(julia_scores) == length(gl_scores) || error("score length mismatch for mode=$mode")
    gl_metrics["gl_max_iter"] <= config.maxiter || error(
        "tutorial mode=$(mode) needed $(Int(round(gl_metrics["gl_max_iter"]))) iterations, above maxiter=$(config.maxiter)"
    )

    abs_diff = abs.(julia_scores .- gl_scores)
    max_abs_diff = isempty(abs_diff) ? 0.0 : maximum(abs_diff)
    return (
        mode = mode,
        n = size(inputs.W, 1),
        nnz = nnz(inputs.W),
        nindustries = length(inputs.industry_labels),
        scenarios = length(gl_scores),
        julia_solve_s = julia_solve_s,
        gl_solve_s = gl_metrics["gl_total_s"],
        speedup_x = gl_metrics["gl_total_s"] / julia_solve_s,
        max_abs_diff = max_abs_diff,
        mean_abs_diff = isempty(abs_diff) ? 0.0 : sum(abs_diff) / length(abs_diff),
        match = max_abs_diff <= config.diff_tol,
        julia_score_range = extrema(julia_scores),
        gl_score_range = (gl_metrics["gl_score_min"], gl_metrics["gl_score_max"]),
        gl_max_iter = Int(round(gl_metrics["gl_max_iter"])),
        julia_threads = Threads.nthreads(),
        julia_threaded = config.julia_threaded,
        gl_ncores = gl_metrics["gl_ncores"] == 0 ? 0 : Int(round(gl_metrics["gl_ncores"])),
        gl_load_balance = gl_metrics["gl_load_balance"] > 0.5,
    )
end

function print_markdown_table(results, config::TutorialBenchmarkConfig)
    gl_mode = config.gl_ncores === false ? "single-core" : "$(config.gl_ncores)-core"
    julia_mode = config.julia_threaded ? "threaded" : "single-core"
    println(
        "solver_tol=$(config.tol) diff_tol=$(config.diff_tol) scenarios=$(config.scenarios) maxiter=$(config.maxiter) julia_mode=$(julia_mode) julia_threads=$(Threads.nthreads()) gl_mode=$(gl_mode) gl_load_balance=$(config.gl_load_balance)"
    )
    println("| mode | firms | nnz | industries | scenarios | julia_solve_s | tutorial_solve_s | julia_speedup_x | max_abs_diff | mean_abs_diff | match<=diff_tol | tutorial_max_iter | julia_threads | gl_ncores |")
    println("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: |")
    for result in results
        println(
            "| ",
            result.mode,
            " | ",
            result.n,
            " | ",
            result.nnz,
            " | ",
            result.nindustries,
            " | ",
            result.scenarios,
            " | ",
            round(result.julia_solve_s; digits = 4),
            " | ",
            round(result.gl_solve_s; digits = 4),
            " | ",
            round(result.speedup_x; digits = 2),
            " | ",
            round(result.max_abs_diff; digits = 6),
            " | ",
            round(result.mean_abs_diff; digits = 6),
            " | ",
            result.match,
            " | ",
            result.gl_max_iter,
            " | ",
            result.julia_threads,
            " | ",
            result.gl_ncores,
            " |",
        )
    end
end

function main(args = ARGS)
    config = parse_args(args)
    inputs = prepare_inputs(_tutorial_root())
    results = map(mode -> compare_mode(inputs, config, mode), config.modes)
    print_markdown_table(results, config)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
