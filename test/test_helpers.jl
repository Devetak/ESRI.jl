deterministic_fixture() = (
    [
        1.0 1.0 0.0
        0.0 1.0 0.0
        0.0 0.0 1.0
    ],
    IndustryInfo([1, 1, 1], [true]),
)

function sparse_fixture(; seed::Int = 123, n::Int = 20, density::Float64 = 0.08, nindustries::Int = 5)
    Random.seed!(seed)
    W = sprand(n, n, density) + 0.1I
    info = IndustryInfo(rand(1:nindustries, n), trues(nindustries))
    return W, info
end

function sample_powerlaw_degree(
    rng,
    n;
    mean_degree::Int = 7,
    alpha::Float64 = 2.3,
    max_degree::Int = min(n, 128),
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
    max_degree::Int = min(n, 128),
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

function powerlaw_fixture(; seed::Int = 42, n::Int = 60, mean_degree::Int = 7, alpha::Float64 = 2.3, nindustries::Int = 6, max_degree::Int = min(n, 128))
    rng = MersenneTwister(seed)
    W, deg = sparse_powerlaw_network(rng, n; mean_degree = mean_degree, alpha = alpha, max_degree = max_degree)
    info = IndustryInfo(rand(rng, 1:nindustries, n), rand(rng, Bool, nindustries))
    return W, info, deg
end

function reference_lockstep_scenario(
    econ::ESRIEconomy{T},
    firm_idx::Integer;
    maxiter::Int = 100,
    tol::Real = 1e-2,
    combine::Symbol = :min,
    shock::Union{Nothing,AbstractVector{<:Real}} = nothing,
) where {T}
    workspace = ESRI._allocate_workspace(T, econ.n, ESRI.num_industries(econ.info))
    ESRI._prepare_shock!(workspace.psi, firm_idx, shock)
    fill!(workspace.previous_upstream, one(T))
    fill!(workspace.previous_downstream, one(T))

    for _ in 1:maxiter
        copyto!(workspace.current_downstream, workspace.previous_downstream)
        ESRI.downstream_shock!(
            econ.downstream_impact_essential,
            econ.downstream_impact_nonessential,
            econ.info,
            econ.row_sums,
            workspace.psi,
            workspace.sigmas,
            workspace.essential_matrix,
            workspace.nonessential_vector,
            workspace.current_downstream,
            workspace.temp_sums,
        )
        ESRI.upstream_step!(
            workspace.current_upstream,
            econ.upstream_impact,
            workspace.previous_upstream,
            workspace.psi,
            econ.row_sums,
        )

        distance = max(
            ESRI._linf_distance(workspace.current_upstream, workspace.previous_upstream),
            ESRI._linf_distance(workspace.current_downstream, workspace.previous_downstream),
        )
        if distance < tol
            break
        end
        copyto!(workspace.previous_upstream, workspace.current_upstream)
        copyto!(workspace.previous_downstream, workspace.current_downstream)
    end

    value = ESRI._reduce_esri(workspace.current_upstream, workspace.current_downstream, econ.row_sums, combine)
    return (
        esri = ESRI._normalize_esri(value, econ),
        upstream = copy(workspace.current_upstream),
        downstream = copy(workspace.current_downstream),
    )
end

assert_bounded(xs; atol = 1e-6) = @test all(0 .<= xs .<= 1 .+ atol)
