include("compare_economic_systemic_risk.jl")

struct AllSolverConfig
    scenarios::Int
    tol::Float64
    diff_tol::Float64
    maxiter::Int
    modes::Vector{Symbol}
end

function parse_all_args(args)
    scenarios = 500
    tol = 1e-2
    diff_tol = 1e-1
    maxiter = 100
    modes = Symbol[:linear, :leontief]

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
        elseif startswith(arg, "--")
            error("unknown argument: $arg")
        else
            scenarios = parse(Int, arg)
        end
        i += 1
    end

    scenarios > 0 || throw(ArgumentError("scenarios must be positive"))
    tol == 1e-2 || throw(ArgumentError("all three solvers can only be aligned at tol=1e-2"))
    diff_tol > 0 || throw(ArgumentError("diff-tol must be positive"))
    maxiter > 0 || throw(ArgumentError("maxiter must be positive"))
    all(in((:linear, :leontief)), modes) || throw(ArgumentError("mode must be linear or leontief"))

    return AllSolverConfig(scenarios, tol, diff_tol, maxiter, unique(modes))
end

function _mean_abs_diff(a, b)
    d = abs.(a .- b)
    return isempty(d) ? 0.0 : sum(d) / length(d)
end

function compare_all_mode(inputs, config::AllSolverConfig, mode::Symbol)
    julia_solve_s, julia_scores = run_julia_mode(
        inputs.W,
        inputs.industries,
        inputs.row_weights,
        mode,
        config.scenarios,
        config.tol,
        config.maxiter,
    )

    tutorial_metrics, tutorial_scores = run_glcascade(inputs.prepared_dir, mode, config.tol, config.scenarios)
    tutorial_max_abs_diff = maximum(abs.(julia_scores .- tutorial_scores))
    tutorial_mean_abs_diff = _mean_abs_diff(julia_scores, tutorial_scores)

    legacy_cfg = LegacyBenchmarkConfig(config.scenarios, config.tol, config.diff_tol, config.maxiter)
    legacy_inputs = write_legacy_inputs(inputs.prepared_dir, inputs.industries, inputs.industry_labels, config.scenarios, mode)
    legacy = run_legacy_model(legacy_inputs, legacy_cfg)
    legacy_max_abs_diff = maximum(abs.(julia_scores .- legacy.scores))
    legacy_mean_abs_diff = _mean_abs_diff(julia_scores, legacy.scores)

    return (
        mode = mode,
        n = size(inputs.W, 1),
        nnz = nnz(inputs.W),
        nindustries = length(inputs.industry_labels),
        scenarios = length(julia_scores),
        julia_solve_s = julia_solve_s,
        tutorial_solve_s = tutorial_metrics["gl_total_s"],
        legacy_solve_s = legacy.elapsed,
        tutorial_speedup_x = tutorial_metrics["gl_total_s"] / julia_solve_s,
        legacy_speedup_x = legacy.elapsed / julia_solve_s,
        tutorial_max_abs_diff = tutorial_max_abs_diff,
        tutorial_mean_abs_diff = tutorial_mean_abs_diff,
        tutorial_match = tutorial_max_abs_diff <= config.diff_tol,
        tutorial_max_iter = Int(round(tutorial_metrics["gl_max_iter"])),
        tutorial_ncores = tutorial_metrics["gl_ncores"] == 0 ? 0 : Int(round(tutorial_metrics["gl_ncores"])),
        legacy_max_abs_diff = legacy_max_abs_diff,
        legacy_mean_abs_diff = legacy_mean_abs_diff,
        legacy_match = legacy_max_abs_diff <= config.diff_tol,
        legacy_max_iter = legacy.max_iter,
    )
end

function print_all_table(results, config::AllSolverConfig)
    println("solver_tol=$(config.tol) diff_tol=$(config.diff_tol) scenarios=$(config.scenarios) maxiter=$(config.maxiter) julia_threads=1 tutorial_ncores=0 legacy_threads=1")
    println("| mode | firms | nnz | industries | scenarios | esri_jl_s | tutorial_cpp_s | econsysrisk_s | tut_speedup_x | legacy_speedup_x | tut_max_abs_diff | legacy_max_abs_diff | tut_match<=diff_tol | legacy_match<=diff_tol | tut_max_iter | legacy_max_iter | tutorial_ncores |")
    println("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: | ---: | ---: |")
    for result in results
        println(
            "| ", result.mode,
            " | ", result.n,
            " | ", result.nnz,
            " | ", result.nindustries,
            " | ", result.scenarios,
            " | ", round(result.julia_solve_s; digits = 4),
            " | ", round(result.tutorial_solve_s; digits = 4),
            " | ", round(result.legacy_solve_s; digits = 4),
            " | ", round(result.tutorial_speedup_x; digits = 2),
            " | ", round(result.legacy_speedup_x; digits = 2),
            " | ", round(result.tutorial_max_abs_diff; digits = 6),
            " | ", round(result.legacy_max_abs_diff; digits = 6),
            " | ", result.tutorial_match,
            " | ", result.legacy_match,
            " | ", result.tutorial_max_iter,
            " | ", result.legacy_max_iter,
            " | ", result.tutorial_ncores,
            " |",
        )
    end
end

function main(args = ARGS)
    config = parse_all_args(args)
    inputs = prepare_inputs(_tutorial_root())
    results = map(mode -> compare_all_mode(inputs, config, mode), config.modes)
    print_all_table(results, config)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
