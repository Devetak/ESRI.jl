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
    components in (:none, :upstream, :downstream, :both) || throw(
        ArgumentError("components must be one of :none, :upstream, :downstream, :both"),
    )
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
    length(final_weights) == n ||
        throw(DimensionMismatch("final_weights must have length n"))
    weights = Vector{T}(undef, n)
    zeroT = zero(T)
    @inbounds for i in eachindex(weights)
        raw = final_weights[i]
        if !isfinite(raw) || raw < zeroT
            throw(
                DomainError(
                    final_weights[i],
                    "final_weights values must be finite and nonnegative",
                ),
            )
        end
        weights[i] = T(raw)
    end
    return weights
end

function _fill_shock!(psi::AbstractVector{T}, shock::AbstractVector{<:Real}) where {T}
    length(shock) == length(psi) ||
        throw(DimensionMismatch("shock vector must have length n"))
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
    combine::Symbol,
) where {T}
    oneT = one(T)
    acc = zero(T)
    n = length(current_upstream)
    @inbounds for j = 1:n
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
        lengths[i] = rowptr[i+1] - rowptr[i]
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

_supports_degree_permutation(econ::ESRIEconomy) =
    econ.upstream_impact isa SparseMatrixCSC &&
    econ.downstream_impact_essential isa SparseMatrixCSR &&
    econ.downstream_impact_nonessential isa SparseMatrixCSR

function _degree_desc_permutation(econ::ESRIEconomy)
    out_degree =
        _csr_row_lengths(econ.downstream_impact_essential) .+
        _csr_row_lengths(econ.downstream_impact_nonessential)
    in_degree = _row_counts(econ.upstream_impact.rowval, econ.n)
    score = out_degree .+ in_degree
    return sortperm(1:econ.n; by = i -> (-score[i], i))
end

_permute_sparse_csr(matrix::SparseMatrixCSR, perm::AbstractVector{Int}) =
    sparsecsr(SparseMatrixCSC(matrix)[perm, perm])

function _permute_sparse_economy(econ::ESRIEconomy, perm::AbstractVector{Int})
    return ESRIEconomy(
        IndustryInfo(econ.info.industry_of_firm[perm], econ.info.input_classification),
        econ.upstream_impact[perm, perm],
        _permute_sparse_csr(econ.downstream_impact_essential, perm),
        _permute_sparse_csr(econ.downstream_impact_nonessential, perm),
        econ.column_sums[perm],
        econ.row_sums[perm],
        econ.total_output,
        econ.n,
    )
end

function _unpermute_values(
    permuted_values::AbstractVector{T},
    perm::AbstractVector{Int},
) where {T}
    values = similar(permuted_values)
    @inbounds for i in eachindex(perm)
        values[perm[i]] = permuted_values[i]
    end
    return values
end

function _allocate_workspace(::Type{T}, n::Int, num_inds::Int) where {T}
    return _ESRIWorkspace(
        Vector{T}(undef, n),
        Vector{T}(undef, n),
        Vector{T}(undef, n),
        Vector{T}(undef, n),
        Vector{T}(undef, n),
        zeros(T, n, num_inds),
        Int[],
        Vector{T}(undef, num_inds),
        Vector{T}(undef, n),
        Vector{T}(undef, n),
    )
end

function _default_shock!(psi::AbstractVector{T}, firm_idx::Integer) where {T}
    fill!(psi, one(T))
    psi[firm_idx] = zero(T)
    return psi
end

_prepare_shock!(
    psi::AbstractVector{T},
    firm_idx::Integer,
    shock::Union{Nothing,AbstractVector{<:Real}},
) where {T} = shock === nothing ? _default_shock!(psi, firm_idx) : _fill_shock!(psi, shock)

_normalize_esri(value::T, econ::ESRIEconomy{T}) where {T} =
    econ.total_output > zero(T) ? value / econ.total_output : value

_normalize_esri(value::T, final_weights::AbstractVector{T}) where {T} = begin
    denominator = sum(final_weights)
    denominator > zero(T) ? value / denominator : value
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

_resolve_components(details::Bool, components::Symbol) = begin
    _validate_components(components)
    details ? :both : components
end

function _compute_single_esri!(
    econ::ESRIEconomy{T},
    workspace::_ESRIWorkspace{T},
    final_weights::AbstractVector{T};
    combine::Symbol = :min,
    maxiter::Int = 100,
    tol::Real = 1e-2,
    verbose::Bool = false,
) where {T}
    current_upstream = workspace.current_upstream
    previous_upstream = workspace.previous_upstream
    current_downstream = workspace.current_downstream
    previous_downstream = workspace.previous_downstream
    psi = workspace.psi
    oneT = one(T)
    fill!(previous_upstream, oneT)
    fill!(previous_downstream, oneT)

    for iter = 1:maxiter
        copyto!(current_downstream, previous_downstream)
        _downstream_shock!(econ, workspace)
        downstream_distance = _linf_distance(current_downstream, previous_downstream)

        upstream_step!(
            current_upstream,
            econ.upstream_impact,
            previous_upstream,
            psi,
            econ.row_sums,
        )
        upstream_distance = _linf_distance(current_upstream, previous_upstream)

        if max(upstream_distance, downstream_distance) < tol
            break
        end
        if verbose && (iter % 10 == 0)
            @info "joint iteration" iter max(upstream_distance, downstream_distance)
        end
        copyto!(previous_upstream, current_upstream)
        copyto!(previous_downstream, current_downstream)
    end

    return _reduce_esri(current_upstream, current_downstream, final_weights, combine)
end

function _solve_default_firm_esri(
    econ::ESRIEconomy{T},
    final_weights::AbstractVector{T},
    firm_idx::Integer,
    combine::Symbol,
    maxiter::Int,
    tol::Real,
    workspace::_ESRIWorkspace{T} = _allocate_workspace(
        T,
        econ.n,
        num_industries(econ.info),
    ),
) where {T}
    _default_shock!(workspace.psi, firm_idx)
    value = _compute_single_esri!(
        econ,
        workspace,
        final_weights;
        combine = combine,
        maxiter = maxiter,
        tol = tol,
        verbose = false,
    )
    return _normalize_esri(value, final_weights)
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
    isempty(firm_sel) && return values
    use_threads = threads && Threads.nthreads() > 1

    if use_threads
        if verbose
            @warn "Ignoring `verbose=true` because progress UI is disabled in threaded mode."
        end
        tasks = map(firm_sel) do firm_idx
            Threads.@spawn begin
                value = _solve_default_firm_esri(econ, weights, firm_idx, combine, maxiter, tol)
                return firm_idx, value
            end
        end
        for task in tasks
            firm_idx, value = fetch(task)
            values[firm_idx] = value
        end
    else
        iter_range = verbose ? ProgressBar(firm_sel; total = length(firm_sel)) : firm_sel
        workspace = _allocate_workspace(T, econ.n, num_industries(econ.info))

        for firm_idx in iter_range
            values[firm_idx] = _solve_default_firm_esri(
                econ,
                weights,
                firm_idx,
                combine,
                maxiter,
                tol,
                workspace,
            )
        end
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
        workspace,
        final_weights;
        combine = combine,
        maxiter = maxiter,
        tol = tol,
        verbose = verbose,
    )
    return _package_result(
        components,
        _normalize_esri(value, final_weights),
        workspace.current_upstream,
        workspace.current_downstream,
    )
end

"""
    esri(econ::ESRIEconomy; maxiter=100, tol=1e-2, verbose=false, threads=false,
         firm_indices=nothing, final_weights=nothing, combine=:min)

Compute one full-closure scenario per selected firm and return one ESRI value per firm.
Set `firm_indices` to select firms; the remaining entries stay zero.
`threads=true` distributes selected firms across available Julia threads.
`final_weights` replaces the output weights and normalizes scores by its sum;
a zero total weight returns `0`. `combine` selects `min(upstream, downstream)`,
`upstream`, or `downstream`. Each scenario stops when both component changes are
below `tol` or after `maxiter` iterations. Each score uses the final iterate.
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

    if firm_indices === nothing && !threads && _supports_degree_permutation(econ)
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

    return _economywide_esri(
        econ,
        weights,
        firm_sel,
        threads,
        verbose,
        combine,
        maxiter,
        tol,
    )
end

"""
    esri(econ::ESRIEconomy, firm_idx::Integer; maxiter=100, tol=1e-2, verbose=false,
         details=false, components=:none, final_weights=nothing, combine=:min,
         shock=nothing)

Solve one scenario and return a scalar, a named tuple, or `ESRIResult`.
The default scenario closes `firm_idx`. Supply `shock ∈ [0, 1]^N` to replace
the default closure. `details=true` is shorthand for `components=:both`.
The solver stops when both component changes are below `tol` or after `maxiter`
iterations. Component results contain the final iterate. Custom `final_weights`
normalize the score by their sum; a zero total weight returns `0`.
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
    compute_esri(weight_matrix, info::IndustryInfo; maxiter=100, tol=1e-2, verbose=false,
                 threads=false, firm_indices=nothing, final_weights=nothing, combine=:min)

Build `ESRIEconomy(weight_matrix, info)` and dispatch to `esri(econ; ...)`.
"""
compute_esri(
    weight_matrix,
    info::IndustryInfo;
    maxiter::Int = 100,
    tol::Real = 1e-2,
    verbose::Bool = false,
    threads::Bool = false,
    firm_indices::Union{Nothing,AbstractVector{<:Integer}} = nothing,
    final_weights::Union{Nothing,AbstractVector{<:Real}} = nothing,
    combine::Symbol = :min,
) = esri(
    ESRIEconomy(weight_matrix, info);
    maxiter,
    tol,
    verbose,
    threads,
    firm_indices,
    final_weights,
    combine,
)

"""
    esri_shock(econ::ESRIEconomy, shock::AbstractVector; maxiter=100, tol=1e-2,
               verbose=false, details=false, components=:none,
               final_weights=nothing, combine=:min)

Solve one scenario from a capacity-cap vector `shock ∈ [0, 1]^N`.
Return a scalar, a named tuple, or `ESRIResult`. The solver stops when both
component changes are below `tol` or after `maxiter` iterations. Component
results contain the final iterate. Custom `final_weights` normalize the score by
their sum; a zero total weight returns `0`.
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
    compute_esri_shock(weight_matrix, info::IndustryInfo, shock::AbstractVector;
                       maxiter=100, tol=1e-2, verbose=false, details=false,
                       components=:none, final_weights=nothing, combine=:min)

Build `ESRIEconomy(weight_matrix, info)` and dispatch to `esri_shock`.
"""
compute_esri_shock(
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
) = esri_shock(
    ESRIEconomy(weight_matrix, info),
    shock;
    maxiter,
    tol,
    verbose,
    details,
    components,
    final_weights,
    combine,
)

"""
    compute_esri(weight_matrix, info::IndustryInfo, firm_idx::Integer;
                 maxiter=100, tol=1e-2, verbose=false, details=false,
                 components=:none, final_weights=nothing, combine=:min, shock=nothing)

Build `ESRIEconomy(weight_matrix, info)` and dispatch to `esri(econ, firm_idx; ...)`.
"""
compute_esri(
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
) = esri(
    ESRIEconomy(weight_matrix, info),
    firm_idx;
    maxiter,
    tol,
    verbose,
    details,
    components,
    final_weights,
    combine,
    shock,
)
