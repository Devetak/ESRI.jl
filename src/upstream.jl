@inline function _sparse_column_accumulate(
    rows::AbstractVector{Int},
    vals::AbstractVector{T},
    previous_upstream::AbstractVector{T},
    start_idx::Int,
    stop_idx::Int,
) where {T}
    acc = zero(T)
    @inbounds for idx = start_idx:stop_idx
        acc += vals[idx] * previous_upstream[rows[idx]]
    end
    return acc
end

function upstream_step!(
    current_upstream::AbstractVector{T},
    upstream_impact::AbstractMatrix{T},
    previous_upstream::AbstractVector{T},
    psi::AbstractVector{T},
    row_sums::AbstractVector{T},
) where {T}
    @inbounds for firm_idx in eachindex(current_upstream)
        if row_sums[firm_idx] == zero(T)
            current_upstream[firm_idx] = min(one(T), psi[firm_idx])
            continue
        end
        acc = zero(T)
        for customer_idx in eachindex(previous_upstream)
            val = upstream_impact[customer_idx, firm_idx]
            if val != zero(T)
                acc += val * previous_upstream[customer_idx]
            end
        end
        current_upstream[firm_idx] = min(acc, psi[firm_idx])
    end
    return current_upstream
end

function upstream_step!(
    current_upstream::AbstractVector{T},
    upstream_impact::SparseMatrixCSC{T},
    previous_upstream::AbstractVector{T},
    psi::AbstractVector{T},
    row_sums::AbstractVector{T},
) where {T}
    rows = upstream_impact.rowval
    colptr = upstream_impact.colptr
    vals = upstream_impact.nzval

    @inbounds for firm_idx in eachindex(current_upstream)
        if row_sums[firm_idx] == zero(T)
            current_upstream[firm_idx] = min(one(T), psi[firm_idx])
            continue
        end
        start_idx = colptr[firm_idx]
        stop_idx = colptr[firm_idx+1] - 1
        acc = _sparse_column_accumulate(rows, vals, previous_upstream, start_idx, stop_idx)
        current_upstream[firm_idx] = min(acc, psi[firm_idx])
    end
    return current_upstream
end
