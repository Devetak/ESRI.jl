include("compare_economic_systemic_risk.jl")
include("sparse_powerlaw_esri.jl")

struct SampledSolverConfig
    sizes::Vector{Int}
    scenarios::Int
    tol::Float64
    diff_tol::Float64
    maxiter::Int
    nindustries::Int
    seed::Int
end

function parse_sampled_args(args)
    sizes = Int[1_000, 2_500, 5_000, 7_500, 10_000]
    scenarios = 500
    tol = 1e-2
    diff_tol = 1e-1
    maxiter = 100
    nindustries = 50
    seed = 42

    positional_sizes = Int[]
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
        elseif arg == "--nindustries"
            i += 1
            nindustries = parse(Int, args[i])
        elseif arg == "--seed"
            i += 1
            seed = parse(Int, args[i])
        elseif startswith(arg, "--")
            error("unknown argument: $arg")
        else
            push!(positional_sizes, parse(Int, arg))
        end
        i += 1
    end

    !isempty(positional_sizes) && (sizes = positional_sizes)
    all(>(0), sizes) || throw(ArgumentError("sizes must be positive"))
    scenarios > 0 || throw(ArgumentError("scenarios must be positive"))
    tol == 1e-2 || throw(ArgumentError("all solver comparisons are aligned at tol=1e-2"))
    diff_tol > 0 || throw(ArgumentError("diff-tol must be positive"))
    maxiter > 0 || throw(ArgumentError("maxiter must be positive"))
    nindustries > 0 || throw(ArgumentError("nindustries must be positive"))

    return SampledSolverConfig(sizes, scenarios, tol, diff_tol, maxiter, nindustries, seed)
end

function synthetic_inputs(n::Int; seed::Int, nindustries::Int, mean_degree::Int = 7, alpha::Float64 = 2.3)
    rng = MersenneTwister(seed)
    W, _ = sparse_powerlaw_network(rng, n; mean_degree = mean_degree, alpha = alpha, max_degree = _default_max_degree(n))
    industries = rand(rng, 1:nindustries, n)
    row_weights = vec(sum(W, dims = 2))
    industry_labels = string.(collect(1:nindustries))

    prepared_dir = mktempdir()
    I, J, V = findnz(W)
    writedlm(joinpath(prepared_dir, "edges.tsv"), hcat(I, J, V), '\t')
    writedlm(joinpath(prepared_dir, "industries.tsv"), industries, '\t')
    writedlm(joinpath(prepared_dir, "industry_labels.tsv"), industry_labels, '\t')
    writedlm(joinpath(prepared_dir, "row_weights.tsv"), row_weights, '\t')
    writedlm(
        joinpath(prepared_dir, "metadata.tsv"),
        ["n" string(n); "nnz" string(nnz(W)); "nindustries" string(nindustries)],
        '\t',
    )

    return (
        prepared_dir = prepared_dir,
        metadata = Dict("n" => string(n), "nnz" => string(nnz(W)), "nindustries" => string(nindustries)),
        W = W,
        industries = industries,
        row_weights = row_weights,
        industry_labels = industry_labels,
    )
end

function compare_sampled_mode(n::Int, mode::Symbol, config::SampledSolverConfig)
    inputs = synthetic_inputs(n; seed = config.seed + n + (mode == :leontief ? 10_000 : 0), nindustries = config.nindustries)

    julia_solve_s, julia_scores = run_julia_mode(
        inputs.W,
        inputs.industries,
        inputs.row_weights,
        mode,
        min(config.scenarios, n),
        config.tol,
        config.maxiter,
    )

    tutorial_metrics, tutorial_scores = run_glcascade(inputs.prepared_dir, mode, config.tol, min(config.scenarios, n))
    tutorial_diff = maximum(abs.(julia_scores .- tutorial_scores))

    legacy_cfg = LegacyBenchmarkConfig(min(config.scenarios, n), config.tol, config.diff_tol, config.maxiter)
    legacy_inputs = write_legacy_inputs(inputs.prepared_dir, inputs.industries, inputs.industry_labels, min(config.scenarios, n), mode)
    legacy = run_legacy_model(legacy_inputs, legacy_cfg)

    return (
        type = "$(mode)_$(n)",
        esricascade_time = julia_solve_s,
        vito_time = legacy.elapsed,
        tutorial_time = tutorial_metrics["gl_total_s"],
        maxdiff = tutorial_diff,
    )
end

function print_sampled_table(results)
    println("| Type | ESRIcascade_TIME | VITO_TIME | TUTORIAL_TIME | MAXDIFF |")
    println("| --- | ---: | ---: | ---: | ---: |")
    for result in results
        println(
            "| ",
            result.type,
            " | ",
            round(result.esricascade_time; digits = 4),
            " | ",
            round(result.vito_time; digits = 4),
            " | ",
            round(result.tutorial_time; digits = 4),
            " | ",
            round(result.maxdiff; digits = 6),
            " |",
        )
    end
end

function main(args = ARGS)
    config = parse_sampled_args(args)
    results = NamedTuple[]
    for n in config.sizes
        push!(results, compare_sampled_mode(n, :linear, config))
        push!(results, compare_sampled_mode(n, :leontief, config))
    end
    print_sampled_table(results)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
