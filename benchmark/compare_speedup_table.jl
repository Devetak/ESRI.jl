using DelimitedFiles
using Printf
using Random
using SparseArrays

include("compare_esri_tutorial.jl")
include("sparse_powerlaw_esri.jl")

struct SpeedupTableConfig
    tol::Float64
    maxiter::Int
    probe_scenarios::Int
    full_candidates_per_thread::Int
    synthetic_10k::Int
    synthetic_50k::Int
    mean_degree::Int
    alpha::Float64
    nindustries::Int
    seed::Int
    batch_sizes::Vector{Union{String, Int}}
    table_out::String
    details_out::String
    attempts_out::String
end

function parse_speedup_args(args)
    tol = 1e-2
    maxiter = 20
    probe_scenarios = 100
    full_candidates_per_thread = 3
    synthetic_10k = 10_000
    synthetic_50k = 50_000
    mean_degree = 7
    alpha = 2.3
    nindustries = 50
    seed = 42
    batch_sizes = Union{String, Int}["all", 256, 512, 1024, 2048, 4096]
    table_out = joinpath(@__DIR__, "speedup_table.md")
    details_out = joinpath(@__DIR__, "speedup_table_details.md")
    attempts_out = joinpath(@__DIR__, "speedup_table_attempts.tsv")

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--tol"
            i += 1
            tol = parse(Float64, args[i])
        elseif arg == "--maxiter"
            i += 1
            maxiter = parse(Int, args[i])
        elseif arg == "--probe-scenarios"
            i += 1
            probe_scenarios = parse(Int, args[i])
        elseif arg == "--full-candidates-per-thread"
            i += 1
            full_candidates_per_thread = parse(Int, args[i])
        elseif arg == "--synthetic-10k"
            i += 1
            synthetic_10k = parse(Int, args[i])
        elseif arg == "--synthetic-50k"
            i += 1
            synthetic_50k = parse(Int, args[i])
        elseif arg == "--mean-degree"
            i += 1
            mean_degree = parse(Int, args[i])
        elseif arg == "--alpha"
            i += 1
            alpha = parse(Float64, args[i])
        elseif arg == "--nindustries"
            i += 1
            nindustries = parse(Int, args[i])
        elseif arg == "--seed"
            i += 1
            seed = parse(Int, args[i])
        elseif arg == "--batch-sizes"
            i += 1
            batch_sizes = map(split(args[i], ',')) do token
                lower = lowercase(strip(token))
                lower == "all" ? "all" : parse(Int, lower)
            end
        elseif arg == "--table-out"
            i += 1
            table_out = args[i]
        elseif arg == "--details-out"
            i += 1
            details_out = args[i]
        elseif arg == "--attempts-out"
            i += 1
            attempts_out = args[i]
        else
            error("unknown argument: $arg")
        end
        i += 1
    end

    tol > 0 || throw(ArgumentError("tol must be positive"))
    maxiter > 0 || throw(ArgumentError("maxiter must be positive"))
    probe_scenarios > 0 || throw(ArgumentError("probe-scenarios must be positive"))
    full_candidates_per_thread > 0 || throw(ArgumentError("full-candidates-per-thread must be positive"))
    synthetic_10k > 0 || throw(ArgumentError("synthetic-10k must be positive"))
    synthetic_50k > 0 || throw(ArgumentError("synthetic-50k must be positive"))
    mean_degree > 0 || throw(ArgumentError("mean-degree must be positive"))
    alpha > 1 || throw(ArgumentError("alpha must be greater than 1"))
    nindustries > 0 || throw(ArgumentError("nindustries must be positive"))

    return SpeedupTableConfig(
        tol,
        maxiter,
        probe_scenarios,
        full_candidates_per_thread,
        synthetic_10k,
        synthetic_50k,
        mean_degree,
        alpha,
        nindustries,
        seed,
        batch_sizes,
        table_out,
        details_out,
        attempts_out,
    )
end

function _read_metrics(path)
    metrics = Dict{String, Float64}()
    for row in eachrow(readdlm(path, '\t', Any))
        metrics[string(row[1])] = row[2] isa Number ? Float64(row[2]) : parse(Float64, string(row[2]))
    end
    return metrics
end

function _read_scores(path)
    return Float64.(vec(readdlm(path, '\t', Float64)))
end

function _write_prepared_dir(prepared_dir::AbstractString, W, industries, industry_labels, row_weights, essential_flags)
    I, J, V = findnz(W)
    writedlm(joinpath(prepared_dir, "edges.tsv"), hcat(I, J, V), '\t')
    writedlm(joinpath(prepared_dir, "industries.tsv"), industries, '\t')
    writedlm(joinpath(prepared_dir, "industry_labels.tsv"), industry_labels, '\t')
    writedlm(joinpath(prepared_dir, "row_weights.tsv"), row_weights, '\t')
    writedlm(joinpath(prepared_dir, "essential.tsv"), Int.(essential_flags), '\t')
end

_synthetic_label(n::Int) = n % 1000 == 0 ? "$(div(n, 1000))k synthetic" : "$(n) synthetic"

function prepare_icio_case()
    inputs = prepare_inputs(_tutorial_root())
    essential_flags = trues(length(inputs.industry_labels))
    writedlm(joinpath(inputs.prepared_dir, "essential.tsv"), Int.(essential_flags), '\t')
    return (
        label = "ICIO",
        prepared_dir = inputs.prepared_dir,
        n = size(inputs.W, 1),
        regime = "leontief",
    )
end

function prepare_synthetic_case(n::Int, config::SpeedupTableConfig)
    rng = MersenneTwister(config.seed + n)
    W, _ = sparse_powerlaw_network(
        rng,
        n;
        mean_degree = config.mean_degree,
        alpha = config.alpha,
        max_degree = _default_max_degree(n),
    )
    used_nindustries = min(config.nindustries, n)
    industries = if n <= used_nindustries
        collect(1:used_nindustries)
    else
        assignments = vcat(collect(1:used_nindustries), rand(rng, 1:used_nindustries, n - used_nindustries))
        shuffle!(rng, assignments)
        assignments
    end
    row_weights = vec(sum(W, dims = 2))
    industry_labels = string.(collect(1:used_nindustries))
    essential_flags = [iseven(i) for i in 1:used_nindustries]
    prepared_dir = mktempdir()
    _write_prepared_dir(prepared_dir, W, industries, industry_labels, row_weights, essential_flags)
    return (
        label = _synthetic_label(n),
        prepared_dir = prepared_dir,
        n = n,
        regime = "mixed-essential",
    )
end

function run_esri_full(prepared_dir::AbstractString, config::SpeedupTableConfig; threads::Int, threaded::Bool)
    metrics_path = joinpath(prepared_dir, "esri_$(threads)t_metrics.tsv")
    scores_path = joinpath(prepared_dir, "esri_$(threads)t_scores.tsv")
    runner = joinpath(@__DIR__, "esri_speed_runner.jl")
    cmd = `julia --project=$(dirname(@__DIR__)) --threads $threads $runner $prepared_dir $(string(config.tol)) $(string(config.maxiter)) $(threaded ? "true" : "false") $metrics_path $scores_path`
    run(cmd)
    return _read_metrics(metrics_path), _read_scores(scores_path)
end

function run_fastcascade_full(
    prepared_dir::AbstractString,
    config::SpeedupTableConfig;
    scenario_count::Int,
    p_market_mode::String,
    ncores::Union{Bool, Int},
    load_balance::Bool,
    batch_size::Union{String, Int},
)
    suffix = replace(
        "p=$(p_market_mode)_n=$(ncores === false ? "false" : string(ncores))_lb=$(load_balance)_b=$(batch_size)",
        r"[^A-Za-z0-9=_-]" => "_",
    )
    metrics_path = joinpath(prepared_dir, "fast_$(suffix)_metrics.tsv")
    scores_path = joinpath(prepared_dir, "fast_$(suffix)_scores.tsv")
    runner = joinpath(@__DIR__, "fastcascade_speed_runner.R")
    ncores_arg = ncores === false ? "false" : string(ncores)
    batch_arg = batch_size isa String ? batch_size : string(batch_size)
    load_balance_arg = load_balance ? "true" : "false"
    cmd = `Rscript $runner $prepared_dir $(string(config.tol)) $(string(scenario_count)) $p_market_mode $ncores_arg $load_balance_arg $batch_arg $metrics_path $scores_path`
    run(cmd)
    return _read_metrics(metrics_path), _read_scores(scores_path)
end

function _diff_metrics(reference_scores::AbstractVector{<:Real}, candidate_scores::AbstractVector{<:Real})
    length(reference_scores) == length(candidate_scores) || error("score length mismatch")
    finite_mask = isfinite.(reference_scores) .& isfinite.(candidate_scores)
    compared = count(finite_mask)
    nonfinite_reference = count(!isfinite, reference_scores)
    nonfinite_candidate = count(!isfinite, candidate_scores)
    if compared == 0
        return (match = false, max_abs_diff = Inf, mean_abs_diff = Inf, compared = 0, nonfinite_reference = nonfinite_reference, nonfinite_candidate = nonfinite_candidate)
    end
    diffs = abs.(reference_scores[finite_mask] .- candidate_scores[finite_mask])
    max_abs_diff = maximum(diffs)
    mean_abs_diff = sum(diffs) / length(diffs)
    match = nonfinite_reference == 0 && nonfinite_candidate == 0 && max_abs_diff <= 1e-2
    return (
        match = match,
        max_abs_diff = max_abs_diff,
        mean_abs_diff = mean_abs_diff,
        compared = compared,
        nonfinite_reference = nonfinite_reference,
        nonfinite_candidate = nonfinite_candidate,
    )
end

function _append_attempt!(
    attempts,
    scenario_label::AbstractString,
    phase::AbstractString,
    thread_label::AbstractString,
    p_market_mode::AbstractString,
    ncores,
    load_balance::Bool,
    batch_size,
    scenario_count::Int,
    solve_s,
    diff_metrics,
    status::AbstractString = "ok",
    error_message::AbstractString = "",
)
    push!(attempts, (
        scenario = String(scenario_label),
        phase = String(phase),
        thread_label = String(thread_label),
        p_market_mode = String(p_market_mode),
        fast_ncores = ncores === false ? 0 : Int(ncores),
        load_balance = load_balance,
        batch_size = batch_size isa String ? String(batch_size) : string(batch_size),
        scenario_count = scenario_count,
        fast_solve_s = solve_s,
        match = diff_metrics.match,
        max_abs_diff = diff_metrics.max_abs_diff,
        mean_abs_diff = diff_metrics.mean_abs_diff,
        compared_scores = diff_metrics.compared,
        nonfinite_reference = diff_metrics.nonfinite_reference,
        nonfinite_fastcascade = diff_metrics.nonfinite_candidate,
        status = String(status),
        error_message = replace(String(error_message), '\t' => ' ', '\n' => ' '),
    ))
end

function _failed_diff_metrics(reference_scores::AbstractVector{<:Real})
    return (
        match = false,
        max_abs_diff = Inf,
        mean_abs_diff = Inf,
        compared = 0,
        nonfinite_reference = count(!isfinite, reference_scores),
        nonfinite_candidate = 0,
    )
end

function matching_p_market_modes(case_info, config::SpeedupTableConfig, reference_scores::AbstractVector{<:Real}, attempts)
    probe_count = min(config.probe_scenarios, length(reference_scores))
    matching_modes = String[]
    for p_market_mode in ("false", "p")
        try
            metrics, scores = run_fastcascade_full(
                case_info.prepared_dir,
                config;
                scenario_count = probe_count,
                p_market_mode = p_market_mode,
                ncores = false,
                load_balance = false,
                batch_size = "all",
            )
            diff_metrics = _diff_metrics(reference_scores[1:probe_count], scores)
            _append_attempt!(
                attempts,
                case_info.label,
                "probe",
                "probe",
                p_market_mode,
                false,
                false,
                "all",
                probe_count,
                metrics["fast_total_s"],
                diff_metrics,
            )
            if diff_metrics.match
                push!(matching_modes, p_market_mode)
            end
        catch err
            _append_attempt!(
                attempts,
                case_info.label,
                "probe",
                "probe",
                p_market_mode,
                false,
                false,
                "all",
                probe_count,
                NaN,
                _failed_diff_metrics(reference_scores[1:probe_count]),
                "error",
                sprint(showerror, err),
            )
        end
    end
    return unique(matching_modes)
end

function candidate_configs(config::SpeedupTableConfig, thread_label::String, p_market_modes::Vector{String})
    configs = NamedTuple[]
    if thread_label == "1T"
        for p_market_mode in p_market_modes, batch_size in config.batch_sizes
            push!(configs, (p_market_mode = p_market_mode, ncores = false, load_balance = false, batch_size = batch_size))
        end
    elseif thread_label == "4T"
        for p_market_mode in p_market_modes, load_balance in (false, true), batch_size in config.batch_sizes
            push!(configs, (p_market_mode = p_market_mode, ncores = 4, load_balance = load_balance, batch_size = batch_size))
        end
    else
        error("unknown thread label: $thread_label")
    end
    return configs
end

function fastest_matching_fastcascade(
    case_info,
    config::SpeedupTableConfig,
    thread_label::String,
    reference_scores::AbstractVector{<:Real},
    attempts,
)
    p_market_modes = matching_p_market_modes(case_info, config, reference_scores, attempts)
    isempty(p_market_modes) && error("No matching fastcascade p_market mode found for $(case_info.label)")

    probe_count = min(config.probe_scenarios, length(reference_scores))
    ranked_candidates = NamedTuple[]
    for candidate in candidate_configs(config, thread_label, p_market_modes)
        try
            metrics, scores = run_fastcascade_full(
                case_info.prepared_dir,
                config;
                scenario_count = probe_count,
                p_market_mode = candidate.p_market_mode,
                ncores = candidate.ncores,
                load_balance = candidate.load_balance,
                batch_size = candidate.batch_size,
            )
            diff_metrics = _diff_metrics(reference_scores[1:probe_count], scores)
            _append_attempt!(
                attempts,
                case_info.label,
                "tune_probe",
                thread_label,
                candidate.p_market_mode,
                candidate.ncores,
                candidate.load_balance,
                candidate.batch_size,
                probe_count,
                metrics["fast_total_s"],
                diff_metrics,
            )
            if diff_metrics.match
                push!(ranked_candidates, (
                    probe_s = metrics["fast_total_s"],
                    p_market_mode = candidate.p_market_mode,
                    ncores = candidate.ncores,
                    load_balance = candidate.load_balance,
                    batch_size = candidate.batch_size,
                ))
            end
        catch err
            _append_attempt!(
                attempts,
                case_info.label,
                "tune_probe",
                thread_label,
                candidate.p_market_mode,
                candidate.ncores,
                candidate.load_balance,
                candidate.batch_size,
                probe_count,
                NaN,
                _failed_diff_metrics(reference_scores[1:probe_count]),
                "error",
                sprint(showerror, err),
            )
        end
    end

    isempty(ranked_candidates) && error("No probe-matching fastcascade configuration found for $(case_info.label) in $thread_label")
    sort!(ranked_candidates, by = x -> x.probe_s)

    best = nothing
    shortlisted_count = min(config.full_candidates_per_thread, length(ranked_candidates))

    function evaluate_full_candidate(candidate)
        try
            metrics, scores = run_fastcascade_full(
                case_info.prepared_dir,
                config;
                scenario_count = length(reference_scores),
                p_market_mode = candidate.p_market_mode,
                ncores = candidate.ncores,
                load_balance = candidate.load_balance,
                batch_size = candidate.batch_size,
            )
            diff_metrics = _diff_metrics(reference_scores, scores)
            _append_attempt!(
                attempts,
                case_info.label,
                "full",
                thread_label,
                candidate.p_market_mode,
                candidate.ncores,
                candidate.load_balance,
                candidate.batch_size,
                length(reference_scores),
                metrics["fast_total_s"],
                diff_metrics,
            )
            if diff_metrics.match
                current = (
                    solve_s = metrics["fast_total_s"],
                    setup_s = get(metrics, "fast_setup_s", NaN),
                    p_market_mode = candidate.p_market_mode,
                    ncores = candidate.ncores,
                    load_balance = candidate.load_balance,
                    batch_size = candidate.batch_size,
                    diff = diff_metrics,
                )
                if best === nothing || current.solve_s < best.solve_s
                    best = current
                end
            end
            return diff_metrics.match
        catch err
            _append_attempt!(
                attempts,
                case_info.label,
                "full",
                thread_label,
                candidate.p_market_mode,
                candidate.ncores,
                candidate.load_balance,
                candidate.batch_size,
                length(reference_scores),
                NaN,
                _failed_diff_metrics(reference_scores),
                "error",
                sprint(showerror, err),
            )
            return false
        end
    end

    for candidate in ranked_candidates[1:shortlisted_count]
        evaluate_full_candidate(candidate)
    end

    if best === nothing && shortlisted_count < length(ranked_candidates)
        for candidate in ranked_candidates[(shortlisted_count + 1):end]
            evaluate_full_candidate(candidate) && break
        end
    end

    best === nothing && error("No matching fastcascade configuration found for $(case_info.label) in $thread_label")
    return best
end

function _render_main_table(rows)
    io = IOBuffer()
    println(io, "| Scenario | Number of Firms | Speedup (1T) | Speedup (4T) |")
    println(io, "| --- | ---: | ---: | ---: |")
    for row in rows
        println(
            io,
            "| ",
            row.scenario,
            " | ",
            row.nfirms,
            " | ",
            @sprintf("%.2fx", row.speedup_1t),
            " | ",
            @sprintf("%.2fx", row.speedup_4t),
            " |",
        )
    end
    return String(take!(io))
end

function _render_details_table(rows, config::SpeedupTableConfig)
    io = IOBuffer()
    println(io, "Strict matched speedup benchmark")
    println(io)
    println(io, "Only configurations with `max_abs_diff <= 1e-2` and fully finite outputs are eligible.")
    println(io, "Rows use full-economy ESRI for every firm. `ICIO` uses aligned Leontief mode. Synthetic rows use mixed essentiality with `essential_industry[i] = iseven(i)`.")
    println(io)
    println(io, "| Scenario | Firms | ESRI 1T s | fastcascade 1T s | Speedup 1T | fastcascade 1T config | ESRI 4T s | fastcascade 4T s | Speedup 4T | fastcascade 4T config | Match Max Abs Diff | Match Mean Abs Diff |")
    println(io, "| --- | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | --- | ---: | ---: |")
    for row in rows
        config_1t = "p_market=$(row.fast_1t.p_market_mode), batch=$(row.fast_1t.batch_size)"
        config_4t = "p_market=$(row.fast_4t.p_market_mode), load_balance=$(row.fast_4t.load_balance), batch=$(row.fast_4t.batch_size)"
        println(
            io,
            "| ", row.scenario,
            " | ", row.nfirms,
            " | ", @sprintf("%.4f", row.esri_1t.solve_s),
            " | ", @sprintf("%.4f", row.fast_1t.solve_s),
            " | ", @sprintf("%.2fx", row.speedup_1t),
            " | ", config_1t,
            " | ", @sprintf("%.4f", row.esri_4t.solve_s),
            " | ", @sprintf("%.4f", row.fast_4t.solve_s),
            " | ", @sprintf("%.2fx", row.speedup_4t),
            " | ", config_4t,
            " | ", @sprintf("%.6g", max(row.fast_1t.diff.max_abs_diff, row.fast_4t.diff.max_abs_diff)),
            " | ", @sprintf("%.6g", max(row.fast_1t.diff.mean_abs_diff, row.fast_4t.diff.mean_abs_diff)),
            " |",
        )
    end
    println(io)
    println(io, "Global settings: `tol=$(config.tol)`, `maxiter=$(config.maxiter)`, `probe_scenarios=$(config.probe_scenarios)`, `full_candidates_per_thread=$(config.full_candidates_per_thread)`.")
    return String(take!(io))
end

function _write_attempts_tsv(path::AbstractString, attempts)
    open(path, "w") do io
        println(io, "scenario\tphase\tthread_label\tp_market_mode\tfast_ncores\tload_balance\tbatch_size\tscenario_count\tfast_solve_s\tmatch\tmax_abs_diff\tmean_abs_diff\tcompared_scores\tnonfinite_reference\tnonfinite_fastcascade\tstatus\terror_message")
        for row in attempts
            println(
                io,
                row.scenario, '\t',
                row.phase, '\t',
                row.thread_label, '\t',
                row.p_market_mode, '\t',
                row.fast_ncores, '\t',
                row.load_balance ? 1 : 0, '\t',
                row.batch_size, '\t',
                row.scenario_count, '\t',
                row.fast_solve_s, '\t',
                row.match ? 1 : 0, '\t',
                row.max_abs_diff, '\t',
                row.mean_abs_diff, '\t',
                row.compared_scores, '\t',
                row.nonfinite_reference, '\t',
                row.nonfinite_fastcascade, '\t',
                row.status, '\t',
                row.error_message,
            )
        end
    end
end

function benchmark_case(case_info, config::SpeedupTableConfig, attempts)
    esri_1t_metrics, esri_1t_scores = run_esri_full(case_info.prepared_dir, config; threads = 1, threaded = false)
    esri_4t_metrics, esri_4t_scores = run_esri_full(case_info.prepared_dir, config; threads = 4, threaded = true)
    esri_thread_diff = _diff_metrics(esri_1t_scores, esri_4t_scores)
    esri_thread_diff.match || error("ESRIcascade 1T/4T outputs do not match for $(case_info.label)")

    fast_1t = fastest_matching_fastcascade(case_info, config, "1T", esri_1t_scores, attempts)
    fast_4t = fastest_matching_fastcascade(case_info, config, "4T", esri_1t_scores, attempts)

    return (
        scenario = case_info.label,
        nfirms = case_info.n,
        esri_1t = (solve_s = esri_1t_metrics["esri_solve_s"], build_s = esri_1t_metrics["esri_build_s"]),
        esri_4t = (solve_s = esri_4t_metrics["esri_solve_s"], build_s = esri_4t_metrics["esri_build_s"]),
        fast_1t = fast_1t,
        fast_4t = fast_4t,
        speedup_1t = fast_1t.solve_s / esri_1t_metrics["esri_solve_s"],
        speedup_4t = fast_4t.solve_s / esri_4t_metrics["esri_solve_s"],
    )
end

function main(args = ARGS)
    config = parse_speedup_args(args)
    attempts = NamedTuple[]
    rows = NamedTuple[]
    cases = (
        prepare_icio_case(),
        prepare_synthetic_case(config.synthetic_10k, config),
        prepare_synthetic_case(config.synthetic_50k, config),
    )

    try
        for case_info in cases
            push!(rows, benchmark_case(case_info, config, attempts))
        end
    catch err
        _write_attempts_tsv(config.attempts_out, attempts)
        rethrow(err)
    end

    main_table = _render_main_table(rows)
    details_table = _render_details_table(rows, config)
    open(config.table_out, "w") do io
        write(io, main_table)
    end
    open(config.details_out, "w") do io
        write(io, details_table)
    end
    _write_attempts_tsv(config.attempts_out, attempts)

    println(main_table)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
