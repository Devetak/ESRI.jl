using ESRIcascade
using Random
using SparseArrays
using Statistics

_default_max_degree(n::Int) = min(max(n - 1, 1), 128)

function sample_powerlaw_degree(
    rng,
    n;
    mean_degree::Int = 7,
    alpha::Float64 = 2.3,
    max_degree::Int = _default_max_degree(n),
)
    ks = collect(1:max_degree)
    cdf = cumsum(ks .^ (-alpha))
    cdf ./= last(cdf)

    deg = Vector{Int}(undef, n)
    @inbounds for i in 1:n
        deg[i] = searchsortedfirst(cdf, rand(rng))
    end

    target = n * mean_degree
    total = sum(deg)
    while total < target
        i = rand(rng, 1:n)
        if deg[i] < max_degree
            deg[i] += 1
            total += 1
        end
    end
    while total > target
        i = rand(rng, 1:n)
        if deg[i] > 1
            deg[i] -= 1
            total -= 1
        end
    end
    return deg
end

function sample_unique_customers!(rng, customers::Vector{Int}, seen::Vector{Int}, stamp::Int, k::Int, n::Int)
    seen[stamp] = stamp
    i = 1
    while i <= k
        customer = rand(rng, 1:n)
        if seen[customer] == stamp
            continue
        end
        seen[customer] = stamp
        customers[i] = customer
        i += 1
    end
    return customers
end

function sparse_powerlaw_network(
    rng,
    n;
    mean_degree::Int = 7,
    alpha::Float64 = 2.3,
    max_degree::Int = _default_max_degree(n),
)
    deg = sample_powerlaw_degree(rng, n; mean_degree = mean_degree, alpha = alpha, max_degree = max_degree)
    rows = Int[]
    cols = Int[]
    vals = Float64[]
    customers = Vector{Int}(undef, maximum(deg))
    seen = zeros(Int, n)
    sizehint!(rows, sum(deg))
    sizehint!(cols, sum(deg))
    sizehint!(vals, sum(deg))

    @inbounds for supplier in 1:n
        k = deg[supplier]
        sample_unique_customers!(rng, customers, seen, supplier, k, n)
        for i in 1:k
            push!(rows, supplier)
            push!(cols, customers[i])
            push!(vals, rand(rng))
        end
    end

    return sparse(rows, cols, vals, n, n), deg
end

function p99_degree(deg::AbstractVector{<:Integer})
    sorted = sort(deg)
    return sorted[clamp(ceil(Int, 0.99 * length(sorted)), 1, length(sorted))]
end

function top1pct_edge_share(deg::AbstractVector{<:Integer})
    sorted = sort(deg)
    k = max(1, round(Int, 0.01 * length(sorted)))
    return sum(view(sorted, length(sorted) - k + 1:length(sorted))) / sum(sorted)
end

function generate_benchmark_input(
    n::Int;
    seed::Int = 42,
    mean_degree::Int = 7,
    alpha::Float64 = 2.3,
    nindustries::Int = 50,
    max_degree::Int = _default_max_degree(n),
)
    rng = MersenneTwister(seed)
    W, deg = sparse_powerlaw_network(rng, n; mean_degree = mean_degree, alpha = alpha, max_degree = max_degree)
    info = IndustryInfo(rand(rng, 1:nindustries, n), rand(rng, Bool, nindustries))
    return W, info, deg
end

function benchmark_once(
    n::Int;
    seed::Int = 42,
    mean_degree::Int = 7,
    alpha::Float64 = 2.3,
    nindustries::Int = 50,
    max_degree::Int = _default_max_degree(n),
    maxiter::Int = 30,
    tol::Float64 = 1e-3,
    threaded::Bool = Threads.nthreads() > 1,
    combine::Symbol = :min,
)
    W, info, deg = generate_benchmark_input(
        n;
        seed = seed,
        mean_degree = mean_degree,
        alpha = alpha,
        nindustries = nindustries,
        max_degree = max_degree,
    )
    build_s = @elapsed econ = ESRIEconomy(W, info)
    solve_s = @elapsed scores = esri(econ; maxiter = maxiter, tol = tol, threads = threaded, combine = combine)

    return (
        threads = Threads.nthreads(),
        n = n,
        nnz = nnz(W),
        mean_degree = mean(deg),
        min_degree = minimum(deg),
        max_degree = maximum(deg),
        p99_degree = p99_degree(deg),
        top1pct_edge_share = top1pct_edge_share(deg),
        alpha = alpha,
        combine = combine,
        build_s = build_s,
        solve_s = solve_s,
        scores = scores,
        score_range = extrema(scores),
        W = W,
        info = info,
    )
end

function print_benchmark(result)
    println("threads=", result.threads)
    println("n=", result.n)
    println("nnz=", result.nnz)
    println("mean_degree=", result.mean_degree)
    println("min_degree=", result.min_degree)
    println("max_degree=", result.max_degree)
    println("p99_degree=", result.p99_degree)
    println("top1pct_edge_share=", result.top1pct_edge_share)
    println("alpha=", result.alpha)
    println("build_s=", result.build_s)
    println("solve_s=", result.solve_s)
    println("score_range=", result.score_range)
end

function main(args = ARGS)
    n = isempty(args) ? 50_000 : parse(Int, args[1])
    print_benchmark(benchmark_once(n))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
