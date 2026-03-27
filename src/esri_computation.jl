function _firm_selection(n::Int, firm_indices::Union{Nothing,AbstractVector{<:Integer}})
    firm_sel = firm_indices === nothing ? collect(1:n) : collect(firm_indices)
    if firm_indices !== nothing
        @assert all(1 .<= firm_sel .<= n) "firm_indices must be in 1:n"
        @assert length(unique(firm_sel)) == length(firm_sel) "firm_indices must be unique"
    end
    return firm_sel
end

function _compute_single_esri!(
    econ::ESRIEconomy{T},
    firm_idx::Integer,
    current_upstream::Vector{T},
    previous_upstream::Vector{T},
    visited::Vector{Bool},
    current_downstream::Vector{T},
    previous_downstream::Vector{T},
    sigmas::Vector{T},
    product_matrix::Matrix{T},
    temp_sums::Vector{T},
    final_shock::Vector{T};
    maxiter::Int = 100,
    tol::Real = 1e-2,
    verbose::Bool = false,
) where {T}
    oneT = one(T)
    zeroT = zero(T)

    if iszero(econ.column_sums[firm_idx])
        fill!(current_upstream, oneT)
        current_upstream[firm_idx] = zeroT
    else
        upstream_shock!(
            econ.upstream_impact,
            firm_idx,
            current_upstream,
            previous_upstream,
            visited;
            maxiter = maxiter,
            tol = tol,
            verbose = verbose,
        )
    end

    if iszero(econ.row_sums[firm_idx])
        fill!(current_downstream, oneT)
        current_downstream[firm_idx] = zeroT
    else
        downstream_shock!(
            econ.downstream_impact,
            econ.info,
            econ.row_sums,
            firm_idx,
            sigmas,
            product_matrix,
            current_downstream,
            previous_downstream,
            temp_sums;
            maxiter = maxiter,
            tol = tol,
            verbose = verbose,
        )
    end

    acc = zeroT
    @inbounds for j in 1:econ.n
        f = min(current_upstream[j], current_downstream[j])
        final_shock[j] = f
        acc += econ.row_sums[j] * (oneT - f)
    end
    return acc
end

"""
    esri(econ::ESRIEconomy; maxiter=100, tol=1e-2, verbose=false, threads=false, firm_indices=nothing)

Compute ESRI for all firms or a subset using a precomputed `ESRIEconomy`.
"""
function esri(
    econ::ESRIEconomy{T};
    maxiter::Int = 100,
    tol::Real = 1e-2,
    verbose::Bool = false,
    threads::Bool = false,
    firm_indices::Union{Nothing,AbstractVector{<:Integer}} = nothing,
) where {T}
    n = econ.n
    values = zeros(T, n)
    firm_sel = _firm_selection(n, firm_indices)
    use_threads = threads && Threads.nthreads() > 1
    num_inds = num_industries(econ.info)

    if use_threads
        if verbose
            @warn "Ignoring `verbose=true` because progress UI is disabled in threaded mode."
        end
        # Julia 1.12 can expose thread IDs larger than `nthreads()` (multiple pools).
        # Allocate by `maxthreadid()` so `threadid()` indexing is always safe.
        nt = Threads.maxthreadid()
        current_upstream = [Vector{T}(undef, n) for _ in 1:nt]
        previous_upstream = [Vector{T}(undef, n) for _ in 1:nt]
        visited = [Vector{Bool}(undef, n) for _ in 1:nt]
        current_downstream = [Vector{T}(undef, n) for _ in 1:nt]
        previous_downstream = [Vector{T}(undef, n) for _ in 1:nt]
        sigmas = [Vector{T}(undef, n) for _ in 1:nt]
        product_matrix = [Matrix{T}(undef, n, num_inds) for _ in 1:nt]
        temp_sums = [Vector{T}(undef, num_inds) for _ in 1:nt]
        final_shock = [Vector{T}(undef, n) for _ in 1:nt]

        Threads.@threads for k in eachindex(firm_sel)
            firm_idx = firm_sel[k]
            tid = Threads.threadid()
            values[firm_idx] = _compute_single_esri!(
                econ,
                firm_idx,
                current_upstream[tid],
                previous_upstream[tid],
                visited[tid],
                current_downstream[tid],
                previous_downstream[tid],
                sigmas[tid],
                product_matrix[tid],
                temp_sums[tid],
                final_shock[tid];
                maxiter = maxiter,
                tol = tol,
                verbose = false,
            )
        end
    else
        current_upstream = Vector{T}(undef, n)
        previous_upstream = Vector{T}(undef, n)
        visited = Vector{Bool}(undef, n)
        current_downstream = Vector{T}(undef, n)
        previous_downstream = Vector{T}(undef, n)
        sigmas = Vector{T}(undef, n)
        product_matrix = Matrix{T}(undef, n, num_inds)
        temp_sums = Vector{T}(undef, num_inds)
        final_shock = Vector{T}(undef, n)

        iter_range = verbose ? ProgressBar(firm_sel; total = length(firm_sel)) : firm_sel
        for firm_idx in iter_range
            values[firm_idx] = _compute_single_esri!(
                econ,
                firm_idx,
                current_upstream,
                previous_upstream,
                visited,
                current_downstream,
                previous_downstream,
                sigmas,
                product_matrix,
                temp_sums,
                final_shock;
                maxiter = maxiter,
                tol = tol,
                verbose = false,
            )
        end
    end

    if econ.total_output > zero(T)
        values ./= econ.total_output
    end
    return values
end

"""
    esri(econ::ESRIEconomy, firm_idx::Integer; maxiter=100, tol=1e-2, verbose=false, details=false, components=:none)

Compute ESRI for a single firm while reusing the prepared economy. By default returns
the scalar ESRI value. Set `details=true` or `components=:both` to retrieve both
upstream and downstream vectors.
"""
function esri(
    econ::ESRIEconomy{T},
    firm_idx::Integer;
    maxiter::Int = 100,
    tol::Real = 1e-2,
    verbose::Bool = false,
    details::Bool = false,
    components::Symbol = :none,
) where {T}
    @assert 1 <= firm_idx <= econ.n "firm_idx must be in 1:n"
    @assert components in (:none, :upstream, :downstream, :both) "components must be one of :none, :upstream, :downstream, :both"

    if details
        components = :both
    end

    n = econ.n
    num_inds = num_industries(econ.info)
    current_upstream = Vector{T}(undef, n)
    previous_upstream = Vector{T}(undef, n)
    visited = Vector{Bool}(undef, n)
    current_downstream = Vector{T}(undef, n)
    previous_downstream = Vector{T}(undef, n)
    sigmas = Vector{T}(undef, n)
    product_matrix = Matrix{T}(undef, n, num_inds)
    temp_sums = Vector{T}(undef, num_inds)
    final_shock = Vector{T}(undef, n)

    value = _compute_single_esri!(
        econ,
        firm_idx,
        current_upstream,
        previous_upstream,
        visited,
        current_downstream,
        previous_downstream,
        sigmas,
        product_matrix,
        temp_sums,
        final_shock;
        maxiter = maxiter,
        tol = tol,
        verbose = verbose,
    )
    esri_value = econ.total_output > zero(T) ? value / econ.total_output : value

    if components == :none
        return esri_value
    elseif components == :upstream
        return (esri = esri_value, upstream = copy(current_upstream))
    elseif components == :downstream
        return (esri = esri_value, downstream = copy(current_downstream))
    end

    return ESRIResult(esri_value, copy(current_upstream), copy(current_downstream))
end

"""
    compute_esri(weight_matrix, info::IndustryInfo; maxiter=100, tol=1e-2, verbose=false, threads=false, firm_indices=nothing)

Convenience wrapper around `ESRIEconomy` + `esri`.
"""
function compute_esri(
    weight_matrix,
    info::IndustryInfo;
    maxiter::Int = 100,
    tol::Real = 1e-2,
    verbose::Bool = false,
    threads::Bool = false,
    firm_indices::Union{Nothing,AbstractVector{<:Integer}} = nothing,
)
    econ = ESRIEconomy(weight_matrix, info)
    return esri(
        econ;
        maxiter = maxiter,
        tol = tol,
        verbose = verbose,
        threads = threads,
        firm_indices = firm_indices,
    )
end
