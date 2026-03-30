#!/usr/bin/env julia
using Random
using SparseArrays
using LinearAlgebra: I
using Printf
using DelimitedFiles
using Statistics: mean, median
using Pkg

using ESRI

"""
Simple power-law-ish synthetic directed production network:
- out-degree per firm follows a Pareto-like distribution (heavy tail)
- edge targets are sampled uniformly
- edge weights are random in (0, 1) scaled by the source "size"
"""
function make_powerlaw_sparse_weights(
    n::Int;
    avg_degree::Int = 10,
    pareto_exponent::Real = 2.3,
    self_weight::Real = 0.1,
    seed::Int = 1,
)::Tuple{SparseMatrixCSC{Float64, Int}, Vector{Float64}}
    rng = MersenneTwister(seed)
    T = Float64

    # Sample raw Pareto sizes: k ~ U^( -1/(alpha-1) )
    # This is a simple continuous approximation; we round to degrees afterwards.
    raw = Vector{T}(undef, n)
    inv = 1 / (pareto_exponent - 1)
    @inbounds for i in 1:n
        raw[i] = (rand(rng)^(-inv))
    end
    scale = avg_degree / (mean(raw))

    outdeg = Vector{Int}(undef, n)
    @inbounds for i in 1:n
        d = Int(round(raw[i] * scale))
        outdeg[i] = clamp(d, 0, n - 1)
    end

    # Build COO and convert to CSC.
    # Note: duplicates are allowed; `sparse` will sum them.
    rows = Vector{Int}(undef, sum(outdeg) + n)
    cols = Vector{Int}(undef, sum(outdeg) + n)
    vals = Vector{T}(undef, sum(outdeg) + n)
    pos = 0

    @inbounds for source in 1:n
        d = outdeg[source]
        for _ in 1:d
            target = rand(rng, 1:n)
            pos += 1
            rows[pos] = source
            cols[pos] = target
            vals[pos] = rand(rng) * raw[source]
        end
    end

    # Add self-loops.
    @inbounds for i in 1:n
        pos += 1
        rows[pos] = i
        cols[pos] = i
        vals[pos] = self_weight
    end

    resize!(rows, pos)
    resize!(cols, pos)
    resize!(vals, pos)

    return sparse(rows, cols, vals, n, n), raw
end

function select_stratified_indices(esri_scores::AbstractVector{<:Real}; bottom::Int = 100, top::Int = 100, middle::Int = 800)
    @assert bottom + top + middle <= length(esri_scores)
    perm = sortperm(esri_scores) # ascending
    bottom_idx = perm[1:bottom]
    top_idx = perm[end - top + 1:end]
    middle_idx = perm[bottom + 1 : bottom + middle]
    return unique(vcat(top_idx, bottom_idx, middle_idx))
end

function elapsed_once(f)::Float64
    t = @elapsed f()
    return t
end

function elapsed_once_with_result(f)
    timed = @timed f()
    # Julia may return either a NamedTuple (newer) or tuple-like structure.
    # Normalize to (time, value).
    if timed isa NamedTuple
        return timed.time, timed.value
    end
    return timed[2], timed[1]
end

const _ECONRISK_MODULE = Ref{Any}(nothing)

function get_econrisk_module()
    if _ECONRISK_MODULE[] !== nothing
        return _ECONRISK_MODULE[]
    end
    econ_dir = joinpath(@__DIR__, "..", "Economic-Systemic-Risk")
    @info "Activating Economic-Systemic-Risk env" econ_dir
    orig_proj = Base.active_project()
    Pkg.activate(econ_dir)
    try
        Pkg.instantiate()
    catch err
        # Some branches update Project.toml without a matching Manifest.
        # Resolve once, then instantiate again to make benchmark runs reproducible.
        @warn "Economic-Systemic-Risk instantiate failed; retrying after resolve" exception = (err, catch_backtrace())
        Pkg.resolve()
        Pkg.instantiate()
    end

    mod = Module(:EconomicSystemicRiskBench)
    Base.include(mod, joinpath(econ_dir, "extern.jl"))
    Base.include(mod, joinpath(econ_dir, "functions.jl"))

    _ECONRISK_MODULE[] = mod
    if orig_proj !== nothing
        Pkg.activate(orig_proj)
    end
    return mod
end

function econrisk_esri_with_mode(
    econ_mod,
    esri_fn,
    M,
    A,
    n::Int,
    maxiter::Int,
)
    # Newer Economic-Systemic-Risk branches expose ESRI(M, A, psi_mat, ParsedARGS),
    # while older ones expose ESRI(M, A). Support both for branch-to-branch comparability.
    psi_all = psi_mat_all_firms(n)
    parsed = Dict("tmax" => maxiter, "timeseries" => false)
    try
        df = Base.invokelatest(esri_fn, M, A, psi_all, parsed)
        return extract_econrisk_esri_vector(econ_mod, df, n), "psi_args"
    catch err
        if !(err isa MethodError)
            rethrow(err)
        end
        vec = Base.invokelatest(esri_fn, M, A)
        return Vector{Float64}(vec), "legacy_no_psi_args"
    end
end

function build_econrisk_market(
    econ_mod,
    W,
    industry_ids::AbstractVector{<:Integer},
    essential_industry::AbstractVector{Bool},
)
    n = size(W, 1)
    @assert size(W, 2) == n

    row_sums = vec(sum(W, dims = 2))
    col_sums = vec(sum(W, dims = 1))

    # Economic-Systemic-Risk code is loaded dynamically via `include`, so under Julia 1.12
    # we must use `invokelatest` to avoid world-age issues when calling newly-defined methods.
    market_ctor = Base.invokelatest(getproperty, econ_mod, :Market)
    company_ctor = Base.invokelatest(getproperty, econ_mod, :Company)
    edge_ctor = Base.invokelatest(getproperty, econ_mod, :Edge)

    M = Base.invokelatest(market_ctor)

    # Create companies with pre-filled baseline outputs/inputs.
    @inbounds for i in 1:n
        M.Companies[i] = Base.invokelatest(
            company_ctor,
            string(i),
            i,
            Int(industry_ids[i]),
            edge_ctor[],
            edge_ctor[],
            Float64(row_sums[i]),
            Float64(col_sums[i]),
        )
    end

    I, J, V = findnz(W)
    @inbounds for k in eachindex(V)
        supplier = I[k]
        customer = J[k]
        suppliernace = Int(industry_ids[supplier])
        customernace = Int(industry_ids[customer])
        edge_type = essential_industry[suppliernace] ? 2 : 1

        e = Base.invokelatest(
            edge_ctor,
            supplier,
            customer,
            suppliernace,
            customernace,
            Float64(V[k]),
            edge_type,
        )

        push!(M.Companies[supplier].customers, e)
        push!(M.Companies[customer].suppliers, e)
        push!(M.Edges, e)
    end

    return M
end

function psi_mat_all_firms(n::Int)
    # Each scenario shocks exactly one firm with shocksize=1 (so ESRI uses psi = 1 - 1 = 0).
    return spdiagm(0 => ones(Float64, n))
end

function psi_mat_one_firm(n::Int, firm_idx::Int)
    return sparse([firm_idx], [1], [1.0], n, 1)
end

function _econrisk_extract_columns(df)
    # Run through `invokelatest` from caller to avoid world-age issues with DataFrames methods.
    return df[!, :index], df[!, :esri]
end

function extract_econrisk_esri_vector(econ_mod, df, n::Int)
    idxs, vals = Base.invokelatest(_econrisk_extract_columns, df)
    out = Vector{Float64}(undef, n)
    @inbounds for k in eachindex(vals)
        out[idxs[k]] = Float64(vals[k])
    end
    return out
end

function diff_metrics(src::AbstractVector{<:Real}, other::AbstractVector{<:Real})
    @assert length(src) == length(other)
    abs_diff = abs.(src .- other)
    return (
        max_abs_diff = maximum(abs_diff),
        mean_abs_diff = mean(abs_diff),
        median_abs_diff = median(abs_diff),
    )
end

function main()
    # Transparent default benchmark:
    # - single-thread only
    # - full ESRI (shock every firm once)
    # - compare this repo vs Economic-Systemic-Risk
    # - report timing + max/mean/median absolute output differences
    n = 10_000
    avg_degree = 10
    num_industries = 50
    essential_industries = 5
    seed = 123
    maxiter = 100
    # External implementation converges using a fixed error threshold of 1e-2.
    tol = 0.01
    pilot_maxiter = 8
    pilot_tol = 1e-1
    threads = false # fixed by request: single thread
    out_name = ""

    for (i, a) in enumerate(ARGS)
        if a == "--n"; n = parse(Int, ARGS[i + 1]) end
        if a == "--avg-degree"; avg_degree = parse(Int, ARGS[i + 1]) end
        if a == "--seed"; seed = parse(Int, ARGS[i + 1]) end
        if a == "--maxiter"; maxiter = parse(Int, ARGS[i + 1]) end
        if a == "--tol"; tol = parse(Float64, ARGS[i + 1]) end
        if a == "--output"; out_name = ARGS[i + 1] end
    end

    if Threads.nthreads() != 1
        error("This benchmark is configured for single-thread execution. Run with JULIA_NUM_THREADS=1.")
    end

    Random.seed!(seed)

    essential_industry = falses(num_industries)
    essential_industry[1:essential_industries] .= true

    # Load external implementation once (activation + instantiate is expensive).
    econ_mod = get_econrisk_module()

    results_dir = joinpath(@__DIR__, "..", "results")
    isdir(results_dir) || mkpath(results_dir)
    default_out = @sprintf("esri_single_thread_n%d_seed%d_avgdeg%d.csv", n, seed, avg_degree)
    out_path = joinpath(results_dir, isempty(out_name) ? default_out : out_name)

    @info "Running transparent benchmark" n avg_degree seed maxiter tol threads
    @info "Generating power-law network" n avg_degree seed
    W, _ = make_powerlaw_sparse_weights(n; avg_degree = avg_degree, seed = seed)

    industry_ids = rand(1:num_industries, n)
    info = IndustryInfo(industry_ids, essential_industry)
    econ = ESRIEconomy(W, info)

    # Warm-up for JIT only (single firm).
    _ = esri(econ; firm_indices = [1], maxiter = pilot_maxiter, tol = pilot_tol, verbose = false, threads = threads)

    @info "Computing full ESRI (src)" n maxiter tol threads
    t_src_s, src_esri = elapsed_once_with_result(() -> esri(econ; maxiter = maxiter, tol = tol, verbose = false, threads = threads))

    M = build_econrisk_market(econ_mod, W, industry_ids, essential_industry)
    build_arrays = Base.invokelatest(getproperty, econ_mod, :buildArrays)
    esri_fn = Base.invokelatest(getproperty, econ_mod, :ESRI)
    A = Base.invokelatest(build_arrays, M)

    # Warm-up for JIT only, using a small pilot market (avoid a second full 10k run).
    pilot_n = min(200, n)
    pilot_W, _ = make_powerlaw_sparse_weights(pilot_n; avg_degree = avg_degree, seed = seed + 1)
    pilot_rng = MersenneTwister(seed + 2)
    pilot_industry_ids = rand(pilot_rng, 1:num_industries, pilot_n)
    pilot_info = IndustryInfo(pilot_industry_ids, essential_industry)
    pilot_M = build_econrisk_market(econ_mod, pilot_W, pilot_industry_ids, essential_industry)
    pilot_A = Base.invokelatest(build_arrays, pilot_M)
    _, econ_mode = econrisk_esri_with_mode(econ_mod, esri_fn, pilot_M, pilot_A, pilot_n, min(pilot_maxiter, maxiter))
    @info "Economic-Systemic-Risk API mode" econ_mode

    @info "Computing full ESRI (Economic-Systemic-Risk)" n maxiter econ_mode
    t_econrisk_s, econ_esri = elapsed_once_with_result(() -> first(econrisk_esri_with_mode(econ_mod, esri_fn, M, A, n, maxiter)))

    metrics = diff_metrics(src_esri, econ_esri)

    open(out_path, "w") do io
        println(io, "n,threads,avg_degree,num_industries,essential_industries,maxiter,tol,seed,t_src_s,t_econrisk_s,max_abs_diff,mean_abs_diff,median_abs_diff")
        @printf(
            io,
            "%d,%s,%d,%d,%d,%d,%g,%d,%.6f,%.6f,%.6e,%.6e,%.6e\n",
            n,
            string(threads),
            avg_degree,
            num_industries,
            essential_industries,
            maxiter,
            tol,
            seed,
            t_src_s,
            t_econrisk_s,
            metrics.max_abs_diff,
            metrics.mean_abs_diff,
            metrics.median_abs_diff,
        )
    end

    @info "Benchmark complete" out_path
end

main()

