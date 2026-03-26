#!/usr/bin/env julia
using Random
using SparseArrays
using LinearAlgebra: I
using Printf
using DelimitedFiles
using Statistics: mean

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

function main()
    # Args (simple parsing; defaults are geared toward local/offline benchmarking)
    n = 100_000
    avg_degree = 10
    num_industries = 50
    essential_industries = 5
    seed = 123
    maxiter = 50
    tol = 1e-3
    pilot_maxiter = 8
    pilot_tol = 1e-1
    threads = false
    rank_by = "size" # "size" (fast) or "esri" (slow pilot ranking)
    bottom = 100
    top = 100
    middle = 800
    random_count = 1000

    for (i, a) in enumerate(ARGS)
        if a == "--n"; n = parse(Int, ARGS[i + 1]) end
        if a == "--avg-degree"; avg_degree = parse(Int, ARGS[i + 1]) end
        if a == "--seed"; seed = parse(Int, ARGS[i + 1]) end
        if a == "--maxiter"; maxiter = parse(Int, ARGS[i + 1]) end
        if a == "--tol"; tol = parse(Float64, ARGS[i + 1]) end
        if a == "--threads"; threads = parse(Bool, ARGS[i + 1]) end
        if a == "--rank-by"; rank_by = ARGS[i + 1] end
        if a == "--bottom"; bottom = parse(Int, ARGS[i + 1]) end
        if a == "--top"; top = parse(Int, ARGS[i + 1]) end
        if a == "--middle"; middle = parse(Int, ARGS[i + 1]) end
        if a == "--random-count"; random_count = parse(Int, ARGS[i + 1]) end
    end

    Random.seed!(seed)

    essential_industry = falses(num_industries)
    essential_industry[1:essential_industries] .= true

    @info "Generating power-law network" n avg_degree seed
    W, firm_sizes = make_powerlaw_sparse_weights(n; avg_degree = avg_degree, seed = seed)

    industry_ids = rand(1:num_industries, n)
    info = IndustryInfo(industry_ids, essential_industry)

    # Rank firms for stratified sampling.
    # - Default (`rank_by=\"size\"`): fast proxy based on the synthetic firm sizes used to generate the network.
    # - Optional (`rank_by=\"esri\"`): slow pilot `compute_esri` run on all firms (may be infeasible for very large n).
    if rank_by == "esri"
        @info "Pilot ESRI run for firm ranking" pilot_maxiter pilot_tol
        esri_scores = compute_esri(W, info; maxiter = pilot_maxiter, tol = pilot_tol, verbose = false, threads = threads)
    else
        esri_scores = firm_sizes
    end

    strat_idx = select_stratified_indices(esri_scores; bottom = bottom, top = top, middle = middle)
    strat_size = bottom + top + middle
    @assert length(strat_idx) == strat_size

    # Sample 1000 additional random firms disjoint from stratified set.
    rng = MersenneTwister(seed + 1)
    in_strat = falses(n)
    in_strat[strat_idx] .= true
    in_sample = copy(in_strat)
    random_idx = Vector{Int}(undef, 0)
    while length(random_idx) < random_count
        cand = rand(rng, 1:n)
        if !in_sample[cand]
            push!(random_idx, cand)
            in_sample[cand] = true
        end
    end

    sample_idx = vcat(strat_idx, random_idx)
    @assert length(sample_idx) == strat_size + random_count

    base_idx = sample_idx[1:1]

    @info "Timing" m=length(sample_idx) n=n
    # Measure time for m=1 and m=2000 to extrapolate full-network runtime.
    @info "Warm-up runs (compilation/JIT)" warmup_m1=true warmup_m=length(sample_idx)
    compute_esri(W, info; firm_indices = base_idx, maxiter = maxiter, tol = tol, verbose = false, threads = threads)
    compute_esri(W, info; firm_indices = sample_idx, maxiter = maxiter, tol = tol, verbose = false, threads = threads)

    t_base = elapsed_once(() -> compute_esri(W, info; firm_indices = base_idx, maxiter = maxiter, tol = tol, verbose = false, threads = threads))
    t_m = elapsed_once(() -> compute_esri(W, info; firm_indices = sample_idx, maxiter = maxiter, tol = tol, verbose = false, threads = threads))

    m = length(sample_idx)
    # Linear model: t(m) ~= t(1) + (m-1)*c  => c = (t_m - t_base)/(m-1)
    c = (t_m - t_base) / (m - 1)
    t_full_est = t_base + c * (n - 1)
    if t_full_est < 0
        @warn "Negative extrapolated runtime; likely measurement noise. Clamping to 0." t_full_est
        t_full_est = 0.0
    end

    results_dir = joinpath(@__DIR__, "..", "results")
    isdir(results_dir) || mkpath(results_dir)
    out_path = joinpath(results_dir, @sprintf("esri_bench_n%d_deg%d_seed%d.csv", n, avg_degree, seed))

    open(out_path, "w") do io
        println(io, "n,avg_degree,num_industries,essential_industries,maxiter,tol,threads,seed,rank_by,m,t_base_s,t_m_s,c_s_per_firm,t_full_est_s")
        @printf(io, "%d,%d,%d,%d,%d,%g,%s,%d,%s,%d,%.6f,%.6f,%.6e,%.6f\n",
            n, avg_degree, num_industries, essential_industries, maxiter, tol, string(threads),
            seed, rank_by, m, t_base, t_m, c, t_full_est)
    end

    @info "Benchmark complete" out_path t_base=t_base t_m=t_m t_full_est=t_full_est
end

main()

