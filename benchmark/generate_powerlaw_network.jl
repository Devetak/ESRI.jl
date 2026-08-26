using Printf
using Random

function main(args)
    length(args) == 2 || error("usage: generate_powerlaw_network.jl FIRMS OUTPUT.csv")

    n = parse(Int, args[1])
    output_path = args[2]
    mean_degree = 8
    alpha = 2.3
    max_degree = min(n - 1, 128)
    n > mean_degree || error("FIRMS must exceed $mean_degree")

    rng = MersenneTwister(42)
    degrees = Vector{Int}(undef, n)
    cdf = cumsum(collect(1:max_degree) .^ (-alpha))
    cdf ./= last(cdf)
    for i = 1:n
        degrees[i] = searchsortedfirst(cdf, rand(rng))
    end

    target_edges = n * mean_degree
    total_edges = sum(degrees)
    while total_edges < target_edges
        firm = rand(rng, 1:n)
        if degrees[firm] < max_degree
            degrees[firm] += 1
            total_edges += 1
        end
    end
    while total_edges > target_edges
        firm = rand(rng, 1:n)
        if degrees[firm] > 1
            degrees[firm] -= 1
            total_edges -= 1
        end
    end

    mkpath(dirname(output_path))
    temporary_path = "$output_path.$(getpid()).tmp"
    seen = zeros(Int, n)
    open(temporary_path, "w") do io
        for supplier = 1:n
            seen[supplier] = supplier
            written = 0
            while written < degrees[supplier]
                customer = rand(rng, 1:n)
                seen[customer] == supplier && continue
                seen[customer] = supplier
                @printf(io, "%d,%d,%.17g\n", supplier, customer, 0.1 + rand(rng))
                written += 1
            end
        end
    end
    mv(temporary_path, output_path; force = true)

    @assert sum(degrees) == target_edges
    @assert minimum(degrees) >= 1
    @assert maximum(degrees) <= max_degree
    println(
        "powerlaw_network firms=$n edges=$target_edges mean_out_degree=$mean_degree " *
        "alpha=$alpha max_out_degree=$(maximum(degrees)) seed=42",
    )
end

main(ARGS)
