using DelimitedFiles
using SparseArrays

include("compare_esri_tutorial.jl")

struct LegacyBenchmarkConfig
    scenarios::Int
    tol::Float64
    diff_tol::Float64
    maxiter::Int
end

function parse_legacy_args(args)
    scenarios = 500
    tol = 1e-2
    diff_tol = 1e-2
    maxiter = 100

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
        elseif startswith(arg, "--")
            error("unknown argument: $arg")
        else
            scenarios = parse(Int, arg)
        end
        i += 1
    end

    scenarios > 0 || throw(ArgumentError("scenarios must be positive"))
    maxiter > 0 || throw(ArgumentError("maxiter must be positive"))
    tol == 1e-2 || throw(ArgumentError("Economic-Systemic-Risk hardcodes convergence threshold 1e-2; use --tol 1e-2"))
    diff_tol > 0 || throw(ArgumentError("diff-tol must be positive"))
    return LegacyBenchmarkConfig(scenarios, tol, diff_tol, maxiter)
end

function _legacy_root()
    return joinpath(@__DIR__, "..", "Economic-Systemic-Risk")
end

function _firm_name(i::Int)
    return "firm_$i"
end

function write_legacy_inputs(
    prepared_dir::AbstractString,
    industries::Vector{Int},
    industry_labels::Vector{String},
    scenarios::Int,
    mode::Symbol,
)
    legacy_dir = mktempdir()
    I, J, V = _read_edges(joinpath(prepared_dir, "edges.tsv"))
    n = length(industries)
    k = min(scenarios, n)
    edge_type = mode == :leontief ? 2 : 1

    input_path = joinpath(legacy_dir, "legacy_input.csv")
    open(input_path, "w") do io
        println(io, "supplier,customer,supplierNACE,customerNACE,weight,type")
        # Register isolated firms with zero-weight type-0 self edges so the legacy parser
        # creates a complete company index without affecting propagation.
        for firm in 1:n
            nace = industry_labels[industries[firm]]
            println(io, _firm_name(firm), ",", _firm_name(firm), ",", nace, ",", nace, ",0.0,0")
        end
        @inbounds for idx in eachindex(I)
            supplier = I[idx]
            customer = J[idx]
            weight = V[idx]
            supplier_nace = industry_labels[industries[supplier]]
            customer_nace = industry_labels[industries[customer]]
            println(io, _firm_name(supplier), ",", _firm_name(customer), ",", supplier_nace, ",", customer_nace, ",", weight, ",", edge_type)
        end
    end

    psi_path = joinpath(legacy_dir, "psi.csv")
    open(psi_path, "w") do io
        println(io, "scenario,firm,shocksize")
        for scenario in 1:k
            println(io, scenario, ",", _firm_name(scenario), ",1.0")
        end
    end

    output_path = joinpath(legacy_dir, "legacy_output.csv")
    return (dir = legacy_dir, input = input_path, psi = psi_path, output = output_path, scenarios = k)
end

function run_legacy_model(legacy_inputs, config::LegacyBenchmarkConfig)
    legacy_root = _legacy_root()
    isdir(legacy_root) || error("expected local checkout at $legacy_root")
    runner = joinpath(dirname(@__DIR__), "benchmark", "economic_systemic_risk_runner.jl")
    output_prefix = joinpath(legacy_inputs.dir, "legacy")
    cmd = `julia --project=$legacy_root --threads 1 $runner $legacy_root $(legacy_inputs.input) $(legacy_inputs.psi) $(string(config.maxiter)) $output_prefix`
    run(Cmd(cmd; dir = legacy_inputs.dir))

    metrics = Dict{String, Float64}()
    for row in eachrow(readdlm(output_prefix * "_metrics.tsv", '\t', Any))
        metrics[string(row[1])] = row[2] isa Number ? Float64(row[2]) : parse(Float64, string(row[2]))
    end
    scores = Float64.(_read_vector(output_prefix * "_scores.tsv", Float64))
    return (
        elapsed = metrics["legacy_total_s"],
        scores = scores,
        max_iter = Int(round(metrics["legacy_max_iter"])),
        mean_iter = metrics["legacy_mean_iter"],
    )
end

function compare_legacy(inputs, config::LegacyBenchmarkConfig, mode::Symbol)
    legacy_inputs = write_legacy_inputs(inputs.prepared_dir, inputs.industries, inputs.industry_labels, config.scenarios, mode)
    julia_solve_s, julia_scores = run_julia_mode(
        inputs.W,
        inputs.industries,
        inputs.row_weights,
        mode,
        config.scenarios,
        config.tol,
        config.maxiter,
    )
    legacy = run_legacy_model(legacy_inputs, config)

    length(julia_scores) == length(legacy.scores) || error("score length mismatch")
    legacy.max_iter <= config.maxiter || error("legacy model needed $(legacy.max_iter) iterations, above maxiter=$(config.maxiter)")

    abs_diff = abs.(julia_scores .- legacy.scores)
    max_abs_diff = isempty(abs_diff) ? 0.0 : maximum(abs_diff)
    max_abs_diff <= config.diff_tol || error("score mismatch: max_abs_diff=$max_abs_diff exceeds diff_tol=$(config.diff_tol)")

    return (
        mode = mode,
        n = size(inputs.W, 1),
        nnz = nnz(inputs.W),
        nindustries = length(inputs.industry_labels),
        scenarios = length(legacy.scores),
        julia_solve_s = julia_solve_s,
        legacy_solve_s = legacy.elapsed,
        speedup_x = legacy.elapsed / julia_solve_s,
        max_abs_diff = max_abs_diff,
        mean_abs_diff = isempty(abs_diff) ? 0.0 : sum(abs_diff) / length(abs_diff),
        legacy_max_iter = legacy.max_iter,
    )
end

function print_legacy_markdown(result)
    println("| target | mode | firms | nnz | industries | scenarios | julia_solve_s | legacy_solve_s | julia_speedup_x | max_abs_diff | mean_abs_diff | legacy_max_iter |")
    println("| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    println(
        "| Economic-Systemic-Risk | ",
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
        round(result.legacy_solve_s; digits = 4),
        " | ",
        round(result.speedup_x; digits = 2),
        " | ",
        round(result.max_abs_diff; digits = 6),
        " | ",
        round(result.mean_abs_diff; digits = 6),
        " | ",
        result.legacy_max_iter,
        " |",
    )
end

function main(args = ARGS)
    config = parse_legacy_args(args)
    inputs = prepare_inputs(_tutorial_root())
    for mode in (:linear, :leontief)
        print_legacy_markdown(compare_legacy(inputs, config, mode))
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
