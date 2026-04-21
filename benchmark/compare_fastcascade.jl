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

function _fastcascade_binary_path()
    return joinpath(@__DIR__, "fastcascade_cli")
end

function _ensure_fastcascade_binary()
    binary_path = _fastcascade_binary_path()
    source_path = joinpath(@__DIR__, "fastcascade_cli.cpp")
    if isfile(binary_path) && mtime(binary_path) >= mtime(source_path)
        return binary_path
    end

    include_dir = "/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/RcppArmadillo/include"
    isdir(include_dir) || error("RcppArmadillo headers not found under $include_dir")

    cmd = `clang++ -O3 -std=c++17 -DARMA_USE_CURRENT -I$include_dir $source_path -o $binary_path -framework Accelerate`
    run(cmd)
    return binary_path
end

function _run_diem_compare(edges_path, industries_path, essential_path, tol, output_path)
    binary_path = _ensure_fastcascade_binary()
    run(`$binary_path $edges_path $industries_path $essential_path $(string(tol)) $output_path`)

    metrics = Dict{String, Float64}()
    for row in eachrow(readdlm(output_path, '\t', Any))
        metrics[string(row[1])] = row[2] isa Number ? Float64(row[2]) : parse(Float64, string(row[2]))
    end
    return metrics
end

function compare_once(
    n::Int;
    seed::Int = 42,
    mean_degree::Int = 7,
    alpha::Float64 = 2.3,
    nindustries::Int = 50,
    max_degree::Int = min(n, 128),
    maxiter::Int = 30,
    tol::Float64 = 1e-3,
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
    )

    tempdir_path = mktempdir()
    edges_path, industries_path, essential_path = _write_network_inputs(tempdir_path, result.W, result.info)
    output_path = joinpath(tempdir_path, "diem_metrics.tsv")
    diem_metrics = _run_diem_compare(edges_path, industries_path, essential_path, tol, output_path)

    esri_total_s = result.build_s + result.solve_s
    return (
        n = result.n,
        nnz = result.nnz,
        max_degree = result.max_degree,
        p99_degree = result.p99_degree,
        top1pct_edge_share = result.top1pct_edge_share,
        esri_build_s = result.build_s,
        esri_solve_s = result.solve_s,
        esri_total_s = esri_total_s,
        diem_total_s = diem_metrics["diem_total_s"],
        esri_total_speedup_x = diem_metrics["diem_total_s"] / esri_total_s,
        esri_score_range = result.score_range,
        diem_score_range = (diem_metrics["diem_score_min"], diem_metrics["diem_score_max"]),
    )
end

function print_markdown_table(results)
    println("| firms | nnz | max_degree | p99_degree | top1pct_edge_share | esri_build_s | esri_solve_s | diem_total_s | esri_total_speedup_x |")
    println("| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
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
            round(result.esri_build_s; digits = 4),
            " | ",
            round(result.esri_solve_s; digits = 4),
            " | ",
            round(result.diem_total_s; digits = 4),
            " | ",
            round(result.esri_total_speedup_x; digits = 2),
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
