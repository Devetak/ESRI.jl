function _firm_selection(n::Int, firm_indices::Union{Nothing,AbstractVector{<:Integer}})
    if firm_indices === nothing
        return collect(1:n)
    end

    firm_sel = Vector{Int}(undef, length(firm_indices))
    seen = falses(n)
    @inbounds for i in eachindex(firm_indices)
        idx = Int(firm_indices[i])
        if !(1 <= idx <= n)
            throw(BoundsError(1:n, idx))
        end
        if seen[idx]
            throw(ArgumentError("firm_indices must be unique"))
        end
        seen[idx] = true
        firm_sel[i] = idx
    end

    return firm_sel
end

function _validate_combine(combine::Symbol)
    combine in (:min, :upstream, :downstream) ||
        throw(ArgumentError("combine must be one of :min, :upstream, :downstream"))
    return combine
end

function _validate_components(components::Symbol)
    components in (:none, :upstream, :downstream, :both) ||
        throw(ArgumentError("components must be one of :none, :upstream, :downstream, :both"))
    return components
end

function _coerce_weights(
    final_weights::Union{Nothing,AbstractVector{<:Real}},
    row_sums::AbstractVector{T},
) where {T}
    if final_weights === nothing
        return row_sums
    end
    n = length(row_sums)
    length(final_weights) == n || throw(DimensionMismatch("final_weights must have length n"))
    weights = Vector{T}(undef, n)
    zeroT = zero(T)
    @inbounds for i in eachindex(weights)
        raw = final_weights[i]
        if !isfinite(raw) || raw < zeroT
            throw(DomainError(final_weights[i], "final_weights values must be finite and nonnegative"))
        end
        weights[i] = T(raw)
    end
    return weights
end

function _coerce_shock(shock::AbstractVector{<:Real}, n::Int, ::Type{T}) where {T}
    length(shock) == n || throw(DimensionMismatch("shock vector must have length n"))
    out = Vector{T}(undef, n)
    oneT = one(T)
    zeroT = zero(T)
    @inbounds for i in eachindex(out)
        raw = shock[i]
        if !isfinite(raw) || raw < zeroT || raw > oneT
            throw(DomainError(shock[i], "shock values must be in [0, 1]"))
        end
        out[i] = T(raw)
    end
    return out
end

function _linf_distance(x::AbstractVector{T}, y::AbstractVector{T}) where {T}
    d = zero(T)
    @inbounds for i in eachindex(x, y)
        v = abs(x[i] - y[i])
        if v > d
            d = v
        end
    end
    return d
end

function _reduce_esri(
    current_upstream::AbstractVector{T},
    current_downstream::AbstractVector{T},
    final_weights::AbstractVector{T},
    combine::Symbol
) where {T}
    oneT = one(T)
    acc = zero(T)
    n = length(current_upstream)
    @inbounds for j in 1:n
        f = if combine === :min
            min(current_upstream[j], current_downstream[j])
        elseif combine === :upstream
            current_upstream[j]
        elseif combine === :downstream
            current_downstream[j]
        else
            throw(ArgumentError("combine must be one of :min, :upstream, :downstream"))
        end
        acc += final_weights[j] * (oneT - f)
    end
    return acc
end

function _compute_single_esri!(
    econ::ESRIEconomy{T},
    current_upstream::Vector{T},
    previous_upstream::Vector{T},
    current_downstream::Vector{T},
    previous_downstream::Vector{T},
    sigmas::Vector{T},
    essential_matrix::Matrix{T},
    temp_sums::Vector{T},
    nonessential_vector::Vector{T},
    psi::Vector{T},
    final_weights::AbstractVector{T};
    combine::Symbol = :min,
    maxiter::Int = 100,
    tol::Real = 1e-2,
    verbose::Bool = false,
) where {T}
    oneT = one(T)
    fill!(previous_upstream, oneT)
    fill!(previous_downstream, oneT)

    for iter = 1:maxiter
        copyto!(current_downstream, previous_downstream)
        downstream_shock!(
            econ.downstream_impact_essential,
            econ.downstream_impact_nonessential,
            econ.info,
            econ.row_sums,
            psi,
            sigmas,
            essential_matrix,
            nonessential_vector,
            current_downstream,
            temp_sums,
        )

        upstream_step!(
            current_upstream,
            econ.upstream_impact,
            previous_upstream,
            psi,
            econ.row_sums,
        )

        distance = max(
            _linf_distance(current_upstream, previous_upstream),
            _linf_distance(current_downstream, previous_downstream),
        )
        if distance < tol
            break
        end
        if verbose && (iter % 10 == 0)
            @info "joint iteration" iter distance
        end
        copyto!(previous_upstream, current_upstream)
        copyto!(previous_downstream, current_downstream)
    end

    return _reduce_esri(current_upstream, current_downstream, final_weights, combine)
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
    final_weights::Union{Nothing,AbstractVector{<:Real}} = nothing,
    combine::Symbol = :min,
) where {T}
    _validate_combine(combine)
    n = econ.n
    weights = _coerce_weights(final_weights, econ.row_sums)

    values = zeros(T, n)
    firm_sel = _firm_selection(n, firm_indices)
    use_threads = threads && Threads.nthreads() > 1
    num_inds = num_industries(econ.info)
    oneT = one(T)
    zeroT = zero(T)

    if use_threads
        if verbose
            @warn "Ignoring `verbose=true` because progress UI is disabled in threaded mode."
        end
        nt = Threads.maxthreadid()
        current_upstream = [Vector{T}(undef, n) for _ in 1:nt]
        previous_upstream = [Vector{T}(undef, n) for _ in 1:nt]
        current_downstream = [Vector{T}(undef, n) for _ in 1:nt]
        previous_downstream = [Vector{T}(undef, n) for _ in 1:nt]
        sigmas = [Vector{T}(undef, n) for _ in 1:nt]
        essential_matrix = [Matrix{T}(undef, n, num_inds) for _ in 1:nt]
        temp_sums = [Vector{T}(undef, num_inds) for _ in 1:nt]
        nonessential_vector = [Vector{T}(undef, n) for _ in 1:nt]
        psi = [Vector{T}(undef, n) for _ in 1:nt]

        Threads.@threads for k in eachindex(firm_sel)
            firm_idx = firm_sel[k]
            tid = Threads.threadid()
            # Initialize psi for single-firm shock
            fill!(psi[tid], oneT)
            psi[tid][firm_idx] = zeroT

            values[firm_idx] = _compute_single_esri!(
                econ,
                current_upstream[tid],
                previous_upstream[tid],
                current_downstream[tid],
                previous_downstream[tid],
                sigmas[tid],
                essential_matrix[tid],
                temp_sums[tid],
                nonessential_vector[tid],
                psi[tid],
                weights;
                combine = combine,
                maxiter = maxiter,
                tol = tol,
                verbose = false,
            )
        end
    else
        current_upstream = Vector{T}(undef, n)
        previous_upstream = Vector{T}(undef, n)
        current_downstream = Vector{T}(undef, n)
        previous_downstream = Vector{T}(undef, n)
        sigmas = Vector{T}(undef, n)
        essential_matrix = Matrix{T}(undef, n, num_inds)
        temp_sums = Vector{T}(undef, num_inds)
        nonessential_vector = Vector{T}(undef, n)
        psi = Vector{T}(undef, n)

        iter_range = verbose ? ProgressBar(firm_sel; total = length(firm_sel)) : firm_sel
        for firm_idx in iter_range
            # Initialize psi for single-firm shock
            fill!(psi, oneT)
            psi[firm_idx] = zeroT

            values[firm_idx] = _compute_single_esri!(
                econ,
                current_upstream,
                previous_upstream,
                current_downstream,
                previous_downstream,
                sigmas,
                essential_matrix,
                temp_sums,
                nonessential_vector,
                psi,
                weights;
                combine = combine,
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

Compute ESRI for a single firm. `upstream` / `downstream` components come from the
coupled fixed-point iteration with common shock cap `psi`.
"""
function esri(
    econ::ESRIEconomy{T},
    firm_idx::Integer;
    maxiter::Int = 100,
    tol::Real = 1e-2,
    verbose::Bool = false,
    details::Bool = false,
    components::Symbol = :none,
    final_weights::Union{Nothing,AbstractVector{<:Real}} = nothing,
    combine::Symbol = :min,
    shock::Union{Nothing,AbstractVector{<:Real}} = nothing,
) where {T}
    if !(1 <= firm_idx <= econ.n)
        throw(BoundsError(1:econ.n, firm_idx))
    end
    _validate_components(components)
    _validate_combine(combine)

    if details
        components = :both
    end

    n = econ.n
    weights = _coerce_weights(final_weights, econ.row_sums)

    num_inds = num_industries(econ.info)
    current_upstream = Vector{T}(undef, n)
    previous_upstream = Vector{T}(undef, n)
    current_downstream = Vector{T}(undef, n)
    previous_downstream = Vector{T}(undef, n)
    sigmas = Vector{T}(undef, n)
    essential_matrix = Matrix{T}(undef, n, num_inds)
    temp_sums = Vector{T}(undef, num_inds)
    nonessential_vector = Vector{T}(undef, n)
    psi = Vector{T}(undef, n)

    if shock !== nothing
        copyto!(psi, _coerce_shock(shock, n, T))
    else
        oneT = one(T)
        zeroT = zero(T)
        fill!(psi, oneT)
        psi[firm_idx] = zeroT
    end

    value = _compute_single_esri!(
        econ,
        current_upstream,
        previous_upstream,
        current_downstream,
        previous_downstream,
        sigmas,
        essential_matrix,
        temp_sums,
        nonessential_vector,
        psi,
        weights;
        combine = combine,
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
    final_weights::Union{Nothing,AbstractVector{<:Real}} = nothing,
    combine::Symbol = :min,
)
    econ = ESRIEconomy(weight_matrix, info)
    return esri(
        econ;
        maxiter = maxiter,
        tol = tol,
        verbose = verbose,
        threads = threads,
        firm_indices = firm_indices,
        final_weights = final_weights,
        combine = combine,
    )
end

"""
    esri_shock(econ::ESRIEconomy, shock::AbstractVector; ...)

Compute ESRI for a custom shock vector `shock` (where 0 is a dead firm and 1 is healthy).
"""
function esri_shock(
    econ::ESRIEconomy{T},
    shock::AbstractVector{<:Real};
    maxiter::Int = 100,
    tol::Real = 1e-2,
    verbose::Bool = false,
    details::Bool = false,
    components::Symbol = :none,
    final_weights::Union{Nothing,AbstractVector{<:Real}} = nothing,
    combine::Symbol = :min,
) where {T}
    n = econ.n
    _validate_components(components)
    _validate_combine(combine)

    if details
        components = :both
    end

    weights = _coerce_weights(final_weights, econ.row_sums)

    num_inds = num_industries(econ.info)
    current_upstream = Vector{T}(undef, n)
    previous_upstream = Vector{T}(undef, n)
    current_downstream = Vector{T}(undef, n)
    previous_downstream = Vector{T}(undef, n)
    sigmas = Vector{T}(undef, n)
    essential_matrix = Matrix{T}(undef, n, num_inds)
    temp_sums = Vector{T}(undef, num_inds)
    nonessential_vector = Vector{T}(undef, n)
    psi = Vector{T}(undef, n)

    copyto!(psi, _coerce_shock(shock, n, T))

    value = _compute_single_esri!(
        econ,
        current_upstream,
        previous_upstream,
        current_downstream,
        previous_downstream,
        sigmas,
        essential_matrix,
        temp_sums,
        nonessential_vector,
        psi,
        weights;
        combine = combine,
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
    compute_esri_shock(weight_matrix, info::IndustryInfo, shock::AbstractVector; ...)

Convenience wrapper for custom-shock ESRI computation.
"""
function compute_esri_shock(
    weight_matrix,
    info::IndustryInfo,
    shock::AbstractVector{<:Real};
    maxiter::Int = 100,
    tol::Real = 1e-2,
    verbose::Bool = false,
    details::Bool = false,
    components::Symbol = :none,
    final_weights::Union{Nothing,AbstractVector{<:Real}} = nothing,
    combine::Symbol = :min,
)
    econ = ESRIEconomy(weight_matrix, info)
    return esri_shock(
        econ,
        shock;
        maxiter = maxiter,
        tol = tol,
        verbose = verbose,
        details = details,
        components = components,
        final_weights = final_weights,
        combine = combine,
    )
end

"""
    compute_esri(weight_matrix, info::IndustryInfo, firm_idx::Integer; ...)

Convenience wrapper for single-firm ESRI.
"""
function compute_esri(
    weight_matrix,
    info::IndustryInfo,
    firm_idx::Integer;
    maxiter::Int = 100,
    tol::Real = 1e-2,
    verbose::Bool = false,
    details::Bool = false,
    components::Symbol = :none,
    final_weights::Union{Nothing,AbstractVector{<:Real}} = nothing,
    combine::Symbol = :min,
    shock::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    econ = ESRIEconomy(weight_matrix, info)
    return esri(
        econ,
        firm_idx;
        maxiter = maxiter,
        tol = tol,
        verbose = verbose,
        details = details,
        components = components,
        final_weights = final_weights,
        combine = combine,
        shock = shock,
    )
end
