#!/usr/bin/env julia
#
# Full-network ESRI timing script (no sampling).
#
# Usage:
#   JULIA_NUM_THREADS=12 julia --project=. benchmark/full_esri_timings.jl
#
# This script generates a synthetic power-law sparse production network and runs:
#   compute_esri(W, info; threads=true, maxiter=..., tol=...)
# for N ∈ {50k, 100k, 250k}, reporting minutes and a markdown table.

using Random
using SparseArrays
using Statistics: mean

using ESRI

function make_powerlaw_sparse_weights(
    n::Int;
    avg_degree::Int = 10,
    pareto_exponent::Real = 2.3,
    self_weight::Real = 0.1,
    seed::Int = 123,
)::SparseMatrixCSC{Float64, Int}
    rng = MersenneTwister(seed)

    raw = Vector{Float64}(undef, n)
    inv = 1 / (pareto_exponent - 1)
    @inbounds for i in 1:n
        raw[i] = rand(rng)^(-inv)
    end

    scale = avg_degree / mean(raw)
    outdeg = Vector{Int}(undef, n)
    @inbounds for i in 1:n
        outdeg[i] = clamp(Int(round(raw[i] * scale)), 0, n - 1)
    end

    nn = sum(outdeg) + n
    rows = Vector{Int}(undef, nn)
    cols = Vector{Int}(undef, nn)
    vals = Vector{Float64}(undef, nn)
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

    @inbounds for i in 1:n
        pos += 1
        rows[pos] = i
        cols[pos] = i
        vals[pos] = self_weight
    end

    resize!(rows, pos)
    resize!(cols, pos)
    resize!(vals, pos)

    return sparse(rows, cols, vals, n, n)
end

function build_info(n::Int; num_industries::Int = 50, essential_industries::Int = 5, seed::Int = 123)
    essential_industry = falses(num_industries)
    essential_industry[1:essential_industries] .= true
    rng = MersenneTwister(seed)
    industry_ids = rand(rng, 1:num_industries, n)
    return IndustryInfo(industry_ids, essential_industry)
end

function time_full_esri_minutes(
    n::Int;
    avg_degree::Int,
    seed::Int,
    maxiter::Int,
    tol::Float64,
    num_industries::Int,
    essential_industries::Int,
)
    W = make_powerlaw_sparse_weights(n; avg_degree = avg_degree, seed = seed)
    info = build_info(n; num_industries = num_industries, essential_industries = essential_industries, seed = seed)

    t = @elapsed compute_esri(W, info; maxiter = maxiter, tol = tol, verbose = true, threads = true)
    return t / 60
end

function main()
    Ns = [50_000, 100_000, 250_000]
    avg_degree = 10
    seed = 123
    maxiter = 10
    tol = 1e-2
    num_industries = 50
    essential_industries = 5

    println("Julia threads: ", Threads.nthreads())
    println("Parameters: avg_degree=$(avg_degree) seed=$(seed) maxiter=$(maxiter) tol=$(tol) num_industries=$(num_industries) essential_industries=$(essential_industries)")
    println()

    minutes = Float64[]
    for n in Ns
        println("Running full ESRI for N=$(n)...")
        tmin = time_full_esri_minutes(
            n;
            avg_degree = avg_degree,
            seed = seed,
            maxiter = maxiter,
            tol = tol,
            num_industries = num_industries,
            essential_industries = essential_industries,
        )
        push!(minutes, tmin)
        println("  time_minutes=$(round(tmin; digits=1))")
        println()
    end

    println("### Full ESRI timings (threads=true)")
    println()
    println("| N firms | Julia threads | Time (minutes) |")
    println("|---:|---:|---:|")
    for (n, tmin) in zip(Ns, minutes)
        println("| $(n) | $(Threads.nthreads()) | $(round(tmin; digits=1)) |")
    end
end

main()

