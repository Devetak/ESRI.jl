using Printf
using Random
using Statistics

include("compare_fastcascade.jl")

struct FastcascadeMatrixConfig
    sizes::Vector{Int}
    scenarios::Int
    tols::Vector{Float64}
    maxiter::Int
    nindustries::Int
    mean_degree::Int
    alpha::Float64
    seed::Int
    output_path::Union{Nothing, String}
end

function parse_fastcascade_matrix_args(args)
    sizes = Int[1_000, 10_000, 50_000]
    scenarios = 100
    tols = Float64[1e-2, 1e-3]
    maxiter = 200
    nindustries = 50
    mean_degree = 7
    alpha = 2.3
    seed = 42
    output_path = nothing

    positional_sizes = Int[]
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--tol"
            i += 1
            push!(tols, parse(Float64, args[i]))
        elseif arg == "--scenarios"
            i += 1
            scenarios = parse(Int, args[i])
        elseif arg == "--tols"
            i += 1
            tols = parse.(Float64, split(args[i], ','))
        elseif arg == "--maxiter"
            i += 1
            maxiter = parse(Int, args[i])
        elseif arg == "--nindustries"
            i += 1
            nindustries = parse(Int, args[i])
        elseif arg == "--mean-degree"
            i += 1
            mean_degree = parse(Int, args[i])
        elseif arg == "--alpha"
            i += 1
            alpha = parse(Float64, args[i])
        elseif arg == "--seed"
            i += 1
            seed = parse(Int, args[i])
        elseif arg == "--output"
            i += 1
            output_path = args[i]
        elseif startswith(arg, "--")
            error("unknown argument: $arg")
        else
            push!(positional_sizes, parse(Int, arg))
        end
        i += 1
    end

    !isempty(positional_sizes) && (sizes = positional_sizes)
    sizes = unique(sizes)
    tols = unique(sort(tols; rev = true))

    all(>(0), sizes) || throw(ArgumentError("sizes must be positive"))
    scenarios > 0 || throw(ArgumentError("scenarios must be positive"))
    all(>(0), tols) || throw(ArgumentError("tols must be positive"))
    maxiter > 0 || throw(ArgumentError("maxiter must be positive"))
    nindustries > 0 || throw(ArgumentError("nindustries must be positive"))
    mean_degree > 0 || throw(ArgumentError("mean-degree must be positive"))
    alpha > 1 || throw(ArgumentError("alpha must be greater than 1"))

    return FastcascadeMatrixConfig(sizes, scenarios, tols, maxiter, nindustries, mean_degree, alpha, seed, output_path)
end

function _essential_flags(mode::Symbol, nindustries::Int)
    if mode === :linear
        return falses(nindustries)
    elseif mode === :mixed
        return [iseven(i) for i in 1:nindustries]
    elseif mode === :leontief
        return trues(nindustries)
    end
    throw(ArgumentError("unknown essential regime: $mode"))
end

function _regime_label(mode::Symbol)
    if mode === :linear
        return "All non-essential"
    elseif mode === :mixed
        return "50% essential"
    elseif mode === :leontief
        return "All essential"
    end
    throw(ArgumentError("unknown essential regime: $mode"))
end

function _prepare_case(n::Int, config::FastcascadeMatrixConfig)
    rng = MersenneTwister(config.seed + n)
    W, deg = sparse_powerlaw_network(
        rng,
        n;
        mean_degree = config.mean_degree,
        alpha = config.alpha,
        max_degree = _default_max_degree(n),
    )
    industries = rand(rng, 1:config.nindustries, n)
    return W, deg, industries
end

function _compare_case(
    W,
    deg,
    industries,
    regime::Symbol,
    tol::Float64,
    combine::Symbol,
    config::FastcascadeMatrixConfig,
)
    essential = _essential_flags(regime, config.nindustries)
    info = IndustryInfo(industries, essential)
    econ = ESRIEconomy(W, info)

    tempdir_path = mktempdir()
    edges_path, industries_path, essential_path = _write_network_inputs(tempdir_path, W, info)
    metrics_path = joinpath(tempdir_path, "fastcascade_metrics.tsv")
    scores_path = joinpath(tempdir_path, "fastcascade_scores.tsv")
    fast_metrics, fast_score_matrix = _run_diem_compare(
        edges_path,
        industries_path,
        essential_path,
        tol,
        metrics_path,
        scores_path,
    )

    esri_scores = nothing
    esri_solve_s = @elapsed begin
        esri_scores = esri(
            econ;
            maxiter = config.maxiter,
            tol = tol,
            threads = false,
            combine = combine,
        )
    end

    fast_scores = fast_score_matrix[:, _combine_column_index(combine)]
    abs_diff = abs.(esri_scores .- fast_scores)

    return (
        firms = size(W, 1),
        nnz = nnz(W),
        mean_degree = mean(deg),
        max_degree = maximum(deg),
        industries = config.nindustries,
        essential_regime = _regime_label(regime),
        essential_share = sum(essential) / length(essential),
        tol = tol,
        combine = String(combine),
        esri_solve_s = esri_solve_s,
        fastcascade_solve_s = fast_metrics["diem_total_s"],
        speedup_x = fast_metrics["diem_total_s"] / esri_solve_s,
        compared_scenarios = length(esri_scores),
        max_abs_diff = isempty(abs_diff) ? 0.0 : maximum(abs_diff),
        mean_abs_diff = isempty(abs_diff) ? 0.0 : sum(abs_diff) / length(abs_diff),
    )
end

function run_fastcascade_matrix(config::FastcascadeMatrixConfig)
    regimes = (:linear, :mixed, :leontief)
    combines = (:min, :downstream, :upstream)
    results = NamedTuple[]

    for n in config.sizes
        W, deg, industries = _prepare_case(n, config)
        for regime in regimes
            for tol in config.tols
                fast_metrics = nothing
                fast_score_matrix = nothing

                essential = _essential_flags(regime, config.nindustries)
                info = IndustryInfo(industries, essential)
                econ = ESRIEconomy(W, info)

                tempdir_path = mktempdir()
                edges_path, industries_path, essential_path = _write_network_inputs(tempdir_path, W, info)
                metrics_path = joinpath(tempdir_path, "fastcascade_metrics.tsv")
                scores_path = joinpath(tempdir_path, "fastcascade_scores.tsv")
                fast_metrics, fast_score_matrix = _run_diem_compare(
                    edges_path,
                    industries_path,
                    essential_path,
                    tol,
                    metrics_path,
                    scores_path;
                    scenario_count = min(config.scenarios, size(W, 1)),
                )

                for combine in combines
                    firm_indices = collect(1:min(config.scenarios, size(W, 1)))
                    esri_scores = nothing
                    esri_solve_s = @elapsed begin
                        esri_scores = esri(
                            econ;
                            maxiter = config.maxiter,
                            tol = tol,
                            threads = false,
                            combine = combine,
                            firm_indices = firm_indices,
                        )
                    end

                    esri_selected_scores = esri_scores[firm_indices]
                    fast_scores = fast_score_matrix[:, _combine_column_index(combine)]
                    abs_diff = abs.(esri_selected_scores .- fast_scores)

                    push!(results, (
                        firms = size(W, 1),
                        nnz = nnz(W),
                        mean_degree = mean(deg),
                        max_degree = maximum(deg),
                        industries = config.nindustries,
                        essential_regime = _regime_label(regime),
                        essential_share = sum(essential) / length(essential),
                        tol = tol,
                        combine = String(combine),
                        esri_solve_s = esri_solve_s,
                        fastcascade_solve_s = fast_metrics["diem_total_s"],
                        speedup_x = fast_metrics["diem_total_s"] / esri_solve_s,
                        compared_scenarios = length(firm_indices),
                        max_abs_diff = isempty(abs_diff) ? 0.0 : maximum(abs_diff),
                        mean_abs_diff = isempty(abs_diff) ? 0.0 : sum(abs_diff) / length(abs_diff),
                    ))
                end
            end
        end
    end

    return results
end

function _render_fastcascade_table(results, config::FastcascadeMatrixConfig)
    io = IOBuffer()
    println(io, "Fastcascade vs ESRIcascade supported-subset benchmark")
    println(io)
    println(
        io,
        "Shared sparse power-law networks, single-core runs, one default single-firm shock per firm, default output weights, no substitution, and common reductions `min`, `downstream`, and `upstream`."
    )
    println(
        io,
        "Generator: mean_degree=$(config.mean_degree), alpha=$(config.alpha), industries=$(config.nindustries), scenarios=$(config.scenarios), maxiter=$(config.maxiter), seed=$(config.seed)."
    )
    println(io)
    println(io, "| Firms | NNZ | Mean Deg | Max Deg | Industries | Essential Regime | Essential Share | Tol | Aggregation | ESRIcascade Solve s | fastcascade Solve s | Speedup x | Compared Scenarios | Max Abs Diff | Mean Abs Diff |")
    println(io, "| ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |")
    for row in results
        println(
            io,
            "| ",
            row.firms,
            " | ",
            row.nnz,
            " | ",
            @sprintf("%.2f", row.mean_degree),
            " | ",
            row.max_degree,
            " | ",
            row.industries,
            " | ",
            row.essential_regime,
            " | ",
            @sprintf("%.2f", row.essential_share),
            " | ",
            @sprintf("%.0e", row.tol),
            " | `",
            row.combine,
            "` | ",
            @sprintf("%.4f", row.esri_solve_s),
            " | ",
            @sprintf("%.4f", row.fastcascade_solve_s),
            " | ",
            @sprintf("%.2f", row.speedup_x),
            " | ",
            row.compared_scenarios,
            " | ",
            @sprintf("%.6g", row.max_abs_diff),
            " | ",
            @sprintf("%.6g", row.mean_abs_diff),
            " |",
        )
    end
    return String(take!(io))
end

function main(args = ARGS)
    config = parse_fastcascade_matrix_args(args)
    results = run_fastcascade_matrix(config)
    table = _render_fastcascade_table(results, config)
    println(table)
    if config.output_path !== nothing
        open(config.output_path, "w") do io
            write(io, table)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
