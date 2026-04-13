function _firm_selection(n::Int, firm_indices::Union{Nothing,AbstractVector{<:Integer}})
    if firm_indices === nothing
        return 1:n
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

function _fill_shock!(psi::AbstractVector{T}, shock::AbstractVector{<:Real}) where {T}
    length(shock) == length(psi) || throw(DimensionMismatch("shock vector must have length n"))
    oneT = one(T)
    zeroT = zero(T)
    @inbounds for i in eachindex(psi, shock)
        raw = shock[i]
        if !isfinite(raw) || raw < zeroT || raw > oneT
            throw(DomainError(raw, "shock values must be in [0, 1]"))
        end
        psi[i] = T(raw)
    end
    return psi
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

function _csr_row_lengths(matrix::SparseMatrixCSR)
    rowptr = matrix.rowptr
    lengths = Vector{Int}(undef, length(rowptr) - 1)
    @inbounds for i in eachindex(lengths)
        lengths[i] = rowptr[i + 1] - rowptr[i]
    end
    return lengths
end

function _row_counts(rows::AbstractVector{Int}, n::Int)
    counts = zeros(Int, n)
    @inbounds for row in rows
        counts[row] += 1
    end
    return counts
end

function _supports_degree_permutation(econ::ESRIEconomy)
    return econ.upstream_impact isa SparseMatrixCSC &&
           econ.downstream_impact_essential isa SparseMatrixCSR &&
           econ.downstream_impact_nonessential isa SparseMatrixCSR
end

function _degree_desc_permutation(econ::ESRIEconomy)
    out_degree = _csr_row_lengths(econ.downstream_impact_essential) .+
                 _csr_row_lengths(econ.downstream_impact_nonessential)
    in_degree = _row_counts(econ.upstream_impact.rowval, econ.n)
    score = out_degree .+ in_degree
    return sortperm(1:econ.n; by = i -> (-score[i], i))
end

function _permute_sparse_csr(matrix::SparseMatrixCSR, perm::AbstractVector{Int})
    return sparsecsr(SparseMatrixCSC(matrix)[perm, perm])
end

function _permute_sparse_economy(econ::ESRIEconomy, perm::AbstractVector{Int})
    return ESRIEconomy(
        IndustryInfo(econ.info.industry_of_firm[perm], econ.info.essential_industry),
        econ.upstream_impact[perm, perm],
        _permute_sparse_csr(econ.downstream_impact_essential, perm),
        _permute_sparse_csr(econ.downstream_impact_nonessential, perm),
        econ.column_sums[perm],
        econ.row_sums[perm],
        econ.total_output,
        econ.n,
    )
end

function _unpermute_values(permuted_values::AbstractVector{T}, perm::AbstractVector{Int}) where {T}
    values = similar(permuted_values)
    @inbounds for i in eachindex(perm)
        values[perm[i]] = permuted_values[i]
    end
    return values
end

struct _ESRIWorkspace{T}
    current_upstream::Vector{T}
    previous_upstream::Vector{T}
    current_downstream::Vector{T}
    previous_downstream::Vector{T}
    sigmas::Vector{T}
    essential_matrix::Matrix{T}
    temp_sums::Vector{T}
    nonessential_vector::Vector{T}
    psi::Vector{T}
end

function _allocate_workspace(::Type{T}, n::Int, num_inds::Int) where {T}
    return _ESRIWorkspace(
        Vector{T}(undef, n),
        Vector{T}(undef, n),
        Vector{T}(undef, n),
        Vector{T}(undef, n),
        Vector{T}(undef, n),
        Matrix{T}(undef, n, num_inds),
        Vector{T}(undef, num_inds),
        Vector{T}(undef, n),
        Vector{T}(undef, n),
    )
end

function _thread_workspaces(::Type{T}, n::Int, num_inds::Int, nt::Int) where {T}
    return [_allocate_workspace(T, n, num_inds) for _ in 1:nt]
end

function _default_shock!(psi::AbstractVector{T}, firm_idx::Integer) where {T}
    fill!(psi, one(T))
    psi[firm_idx] = zero(T)
    return psi
end

function _prepare_shock!(
    psi::AbstractVector{T},
    firm_idx::Integer,
    shock::Union{Nothing,AbstractVector{<:Real}},
) where {T}
    return shock === nothing ? _default_shock!(psi, firm_idx) : _fill_shock!(psi, shock)
end

function _normalize_esri(value::T, econ::ESRIEconomy{T}) where {T}
    return econ.total_output > zero(T) ? value / econ.total_output : value
end

function _package_result(
    components::Symbol,
    esri_value,
    current_upstream::AbstractVector,
    current_downstream::AbstractVector,
)
    if components == :none
        return esri_value
    elseif components == :upstream
        return (esri = esri_value, upstream = copy(current_upstream))
    elseif components == :downstream
        return (esri = esri_value, downstream = copy(current_downstream))
    end
    return ESRIResult(esri_value, copy(current_upstream), copy(current_downstream))
end

function _resolve_components(details::Bool, components::Symbol)
    _validate_components(components)
    return details ? :both : components
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
    zeroT = zero(T)
    fill!(previous_upstream, oneT)
    fill!(previous_downstream, oneT)
    downstream_active = true

    for iter = 1:maxiter
        downstream_distance = zeroT
        if downstream_active
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
            downstream_distance = _linf_distance(current_downstream, previous_downstream)
        end

        upstream_step!(
            current_upstream,
            econ.upstream_impact,
            previous_upstream,
            psi,
            econ.row_sums,
        )
        upstream_distance = _linf_distance(current_upstream, previous_upstream)

        if downstream_active
            copyto!(previous_downstream, current_downstream)
            downstream_active = downstream_distance >= tol
        end

        if upstream_distance < tol && !downstream_active
            break
        end
        if verbose && (iter % 10 == 0)
            @info "joint iteration" iter max(upstream_distance, downstream_distance)
        end
        copyto!(previous_upstream, current_upstream)
    end

    return _reduce_esri(current_upstream, current_downstream, final_weights, combine)
end

function _economywide_esri(
    econ::ESRIEconomy{T},
    weights::AbstractVector{T},
    firm_sel,
    threads::Bool,
    verbose::Bool,
    combine::Symbol,
    maxiter::Int,
    tol::Real,
) where {T}
    values = zeros(T, econ.n)
    use_threads = threads && Threads.nthreads() > 1
    num_inds = num_industries(econ.info)

    if use_threads
        if verbose
            @warn "Ignoring `verbose=true` because progress UI is disabled in threaded mode."
        end
        workspaces = _thread_workspaces(T, econ.n, num_inds, Threads.maxthreadid())

        Threads.@threads :static for k in eachindex(firm_sel)
            firm_idx = firm_sel[k]
            workspace = workspaces[Threads.threadid()]
            _default_shock!(workspace.psi, firm_idx)

            values[firm_idx] = _compute_single_esri!(
                econ,
                workspace.current_upstream,
                workspace.previous_upstream,
                workspace.current_downstream,
                workspace.previous_downstream,
                workspace.sigmas,
                workspace.essential_matrix,
                workspace.temp_sums,
                workspace.nonessential_vector,
                workspace.psi,
                weights;
                combine = combine,
                maxiter = maxiter,
                tol = tol,
                verbose = false,
            )
        end
    else
        workspace = _allocate_workspace(T, econ.n, num_inds)
        iter_range = verbose ? ProgressBar(firm_sel; total = length(firm_sel)) : firm_sel

        for firm_idx in iter_range
            _default_shock!(workspace.psi, firm_idx)

            values[firm_idx] = _compute_single_esri!(
                econ,
                workspace.current_upstream,
                workspace.previous_upstream,
                workspace.current_downstream,
                workspace.previous_downstream,
                workspace.sigmas,
                workspace.essential_matrix,
                workspace.temp_sums,
                workspace.nonessential_vector,
                workspace.psi,
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

function _run_scenario!(
    econ::ESRIEconomy{T},
    workspace::_ESRIWorkspace{T},
    final_weights::AbstractVector{T};
    combine::Symbol = :min,
    maxiter::Int = 100,
    tol::Real = 1e-2,
    verbose::Bool = false,
    components::Symbol = :none,
) where {T}
    value = _compute_single_esri!(
        econ,
        workspace.current_upstream,
        workspace.previous_upstream,
        workspace.current_downstream,
        workspace.previous_downstream,
        workspace.sigmas,
        workspace.essential_matrix,
        workspace.temp_sums,
        workspace.nonessential_vector,
        workspace.psi,
        final_weights;
        combine = combine,
        maxiter = maxiter,
        tol = tol,
        verbose = verbose,
    )
    return _package_result(
        components,
        _normalize_esri(value, econ),
        workspace.current_upstream,
        workspace.current_downstream,
    )
end

"""
    esri(econ::ESRIEconomy; maxiter=100, tol=1e-2, verbose=false, threads=false,
         firm_indices=nothing, final_weights=nothing, combine=:min)

Compute one default firm shock per selected firm and return one ESRI value per firm.
Entries outside `firm_indices` stay zero. `final_weights` changes only the numerator,
and `combine` picks `min(upstream, downstream)`, `upstream`, or `downstream`.
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
    firm_sel = _firm_selection(n, firm_indices)

    if firm_indices === nothing && _supports_degree_permutation(econ)
        perm = _degree_desc_permutation(econ)
        permuted_econ = _permute_sparse_economy(econ, perm)
        permuted_values = _economywide_esri(
            permuted_econ,
            weights[perm],
            1:n,
            threads,
            verbose,
            combine,
            maxiter,
            tol,
        )
        return _unpermute_values(permuted_values, perm)
    end

    return _economywide_esri(econ, weights, firm_sel, threads, verbose, combine, maxiter, tol)
end

"""
    esri(econ::ESRIEconomy, firm_idx::Integer; maxiter=100, tol=1e-2, verbose=false,
         details=false, components=:none, final_weights=nothing, combine=:min,
         shock=nothing)

Solve one scenario and return a scalar, a named tuple, or `ESRIResult`.
By default the scenario closes `firm_idx`. If `shock` is given, it must lie in `[0, 1]^N`
and replaces the default closure. `details=true` is shorthand for `components=:both`.
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
    components = _resolve_components(details, components)
    _validate_combine(combine)

    n = econ.n
    weights = _coerce_weights(final_weights, econ.row_sums)

    workspace = _allocate_workspace(T, n, num_industries(econ.info))
    _prepare_shock!(workspace.psi, firm_idx, shock)
    return _run_scenario!(
        econ,
        workspace,
        weights;
        combine = combine,
        maxiter = maxiter,
        tol = tol,
        verbose = verbose,
        components = components,
    )
end

"""
    compute_esri(weight_matrix, info::IndustryInfo; maxiter=100, tol=1e-2, verbose=false, threads=false, firm_indices=nothing)

Build `ESRIEconomy(weight_matrix, info)` and dispatch to `esri(econ; ...)`.
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

Solve one scenario from a capacity-cap vector `shock ∈ [0, 1]^N`.
Return a scalar, a named tuple, or `ESRIResult`.
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
    components = _resolve_components(details, components)
    _validate_combine(combine)

    weights = _coerce_weights(final_weights, econ.row_sums)
    workspace = _allocate_workspace(T, n, num_industries(econ.info))
    _fill_shock!(workspace.psi, shock)
    return _run_scenario!(
        econ,
        workspace,
        weights;
        combine = combine,
        maxiter = maxiter,
        tol = tol,
        verbose = verbose,
        components = components,
    )
end

"""
    compute_esri_shock(weight_matrix, info::IndustryInfo, shock::AbstractVector; ...)

Build `ESRIEconomy(weight_matrix, info)` and dispatch to `esri_shock`.
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

Build `ESRIEconomy(weight_matrix, info)` and dispatch to `esri(econ, firm_idx; ...)`.
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
