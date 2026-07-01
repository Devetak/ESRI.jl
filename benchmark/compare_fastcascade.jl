using DelimitedFiles

include("sparse_powerlaw_esri.jl")

function _write_network_inputs(tempdir_path::AbstractString, W, info)
    edges_path = joinpath(tempdir_path, "edges.tsv")
    industries_path = joinpath(tempdir_path, "industries.tsv")
    essential_path = joinpath(tempdir_path, "essential.tsv")

    I, J, V = findnz(W)
    writedlm(edges_path, hcat(I, J, V), '\t')
    writedlm(industries_path, info.industry_of_firm, '\t')
    writedlm(essential_path, Int.(info.essential_industry), '\t')
    return edges_path, industries_path, essential_path
end

function _fastcascade_runner_path()
    return joinpath(@__DIR__, "fastcascade_runner.R")
end

function _read_metrics(path)
    metrics = Dict{String, Float64}()
    for row in eachrow(readdlm(path, '\t', Any))
        metrics[string(row[1])] = row[2] isa Number ? Float64(row[2]) : parse(Float64, string(row[2]))
    end
    return metrics
end

function _read_score_matrix(path)
    scores = readdlm(path, '\t', Float64)
    if ndims(scores) == 1
        return reshape(Float64.(scores), :, 1)
    end
    return Float64.(scores)
end

function _combine_column_index(combine::Symbol)
    if combine === :min
        return 1
    elseif combine === :downstream
        return 2
    elseif combine === :upstream
        return 3
    end
    throw(ArgumentError("combine must be one of :min, :downstream, :upstream"))
end

function _run_diem_compare(edges_path, industries_path, essential_path, tol, metrics_path, scores_path; scenario_count::Union{Nothing,Int} = nothing)
    runner_path = _fastcascade_runner_path()
    scenario_arg = scenario_count === nothing ? "all" : string(scenario_count)
    run(`Rscript $runner_path $edges_path $industries_path $essential_path $(string(tol)) $scenario_arg $metrics_path $scores_path`)
    return _read_metrics(metrics_path), _read_score_matrix(scores_path)
end

function compare_once(
    n::Int;
    seed::Int = 42,
    mean_degree::Int = 7,
    alpha::Float64 = 2.3,
    nindustries::Int = 50,
    max_degree::Int = _default_max_degree(n),
    maxiter::Int = 10,
    tol::Float64 = 1e-2,
    combine::Symbol = :min,
)
    result = benchmark_once(
        n;
        seed = seed,
        mean_degree = mean_degree,
        alpha = alpha,
        nindustries = nindustries,
        max_degree = max_degree,
        maxiter = maxiter,
        tol = tol,
        threaded = false,
        combine = combine,
    )

    tempdir_path = mktempdir()
    edges_path, industries_path, essential_path = _write_network_inputs(tempdir_path, result.W, result.info)
    metrics_path = joinpath(tempdir_path, "diem_metrics.tsv")
    scores_path = joinpath(tempdir_path, "diem_scores.tsv")
    diem_metrics, diem_score_matrix = _run_diem_compare(
        edges_path,
        industries_path,
        essential_path,
        tol,
        metrics_path,
        scores_path,
    )
    diem_scores = diem_score_matrix[:, _combine_column_index(combine)]

    finite_mask = isfinite.(result.scores) .& isfinite.(diem_scores)
    finite_diffs = abs.(result.scores[finite_mask] .- diem_scores[finite_mask])
    max_abs_diff = isempty(finite_diffs) ? NaN : maximum(finite_diffs)
    mean_abs_diff = isempty(finite_diffs) ? NaN : sum(finite_diffs) / length(finite_diffs)

    return (
        n = result.n,
        nnz = result.nnz,
        max_degree = result.max_degree,
        p99_degree = result.p99_degree,
        top1pct_edge_share = result.top1pct_edge_share,
        combine = combine,
        esri_solve_s = result.solve_s,
        diem_solve_s = diem_metrics["diem_total_s"],
        esri_solve_speedup_x = diem_metrics["diem_total_s"] / result.solve_s,
        compared_scores = count(finite_mask),
        diem_nonfinite_scores = count(!isfinite, diem_scores),
        max_abs_diff = max_abs_diff,
        mean_abs_diff = mean_abs_diff,
        esri_score_range = result.score_range,
        diem_score_range = (diem_metrics["diem_score_min"], diem_metrics["diem_score_max"]),
    )
end

function print_markdown_table(results)
    println("| firms | nnz | max_degree | p99_degree | top1pct_edge_share | esri_solve_s | diem_solve_s | esri_solve_speedup_x | compared_scores | diem_nonfinite_scores | max_abs_diff | mean_abs_diff |")
    println("| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for result in results
        println(
            "| ",
            result.n,
            " | ",
            result.nnz,
            " | ",
            result.max_degree,
            " | ",
            result.p99_degree,
            " | ",
            round(result.top1pct_edge_share; digits = 4),
            " | ",
            round(result.esri_solve_s; digits = 4),
            " | ",
            round(result.diem_solve_s; digits = 4),
            " | ",
            round(result.esri_solve_speedup_x; digits = 2),
            " | ",
            result.compared_scores,
            " | ",
            result.diem_nonfinite_scores,
            " | ",
            round(result.max_abs_diff; digits = 6),
            " | ",
            round(result.mean_abs_diff; digits = 6),
            " |",
        )
    end
end

function main(args = ARGS)
    sizes = isempty(args) ? [1_000, 10_000, 50_000, 100_000] : parse.(Int, args)
    results = map(compare_once, sizes)
    print_markdown_table(results)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
