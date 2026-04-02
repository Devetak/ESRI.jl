"""
    create_upstream_impact_matrix(weight_matrix)

Build the normalized upstream impact matrix.
"""
function create_upstream_impact_matrix(weight_matrix::AbstractMatrix{T}) where {T<:Real}
    nrows, ncols = size(weight_matrix)
    result = zeros(eltype(weight_matrix), nrows, ncols)
    row_sums = vec(sum(weight_matrix, dims = 2))
    @inbounds for source = 1:nrows
        denom = row_sums[source]
        if denom == 0
            continue
        end
        for target = 1:ncols
            val = weight_matrix[source, target]
            if val != 0
                result[target, source] = val / denom
            end
        end
    end
    return result
end

function create_upstream_impact_matrix(weight_matrix::SparseMatrixCSC{T}) where {T<:Real}
    nrows, ncols = size(weight_matrix)
    rows = weight_matrix.rowval
    colptr = weight_matrix.colptr
    vals = weight_matrix.nzval

    row_sums = zeros(T, nrows)
    @inbounds for col = 1:ncols
        for idx = colptr[col]:(colptr[col+1]-1)
            row = rows[idx]
            row_sums[row] += vals[idx]
        end
    end

    nnz = length(vals)
    new_rows = Vector{Int}(undef, nnz)
    new_cols = Vector{Int}(undef, nnz)
    new_vals = Vector{T}(undef, nnz)

    idx_count = 0
    @inbounds for col = 1:ncols
        for idx = colptr[col]:(colptr[col+1]-1)
            source = rows[idx]
            target = col
            val = vals[idx]
            denom = row_sums[source]
            if denom != 0
                idx_count += 1
                new_rows[idx_count] = target
                new_cols[idx_count] = source
                new_vals[idx_count] = val / denom
            end
        end
    end

    resize!(new_rows, idx_count)
    resize!(new_cols, idx_count)
    resize!(new_vals, idx_count)

    return sparse(new_rows, new_cols, new_vals, ncols, nrows)
end

function compute_downstream_impact_matrices(
    weight_matrix::AbstractMatrix{T},
    info::IndustryInfo,
) where {T<:Real}
    nrows, ncols = size(weight_matrix)
    essential = zeros(T, nrows, ncols)
    nonessential = zeros(T, nrows, ncols)
    num_inds = num_industries(info)
    essential_by_industry = zeros(T, num_inds)

    @inbounds for col = 1:ncols
        fill!(essential_by_industry, zero(T))
        all_suppliers_total = zero(T)
        for idx = 1:nrows
            val = weight_matrix[idx, col]
            if val == 0
                continue
            end
            all_suppliers_total += val
            if is_essential(info, idx)
                essential_by_industry[get_industry(info, idx)] += val
            end
        end

        for idx = 1:nrows
            val = weight_matrix[idx, col]
            if val == 0
                continue
            end
            if is_essential(info, idx)
                denom = essential_by_industry[get_industry(info, idx)]
                essential[idx, col] = denom == 0 ? zero(T) : val / denom
            else
                nonessential[idx, col] = all_suppliers_total == 0 ? zero(T) : val / all_suppliers_total
            end
        end
    end

    return essential, nonessential
end

function compute_downstream_impact_matrices(
    weight_matrix::SparseMatrixCSC{T},
    info::IndustryInfo,
) where {T<:Real}
    nrows, ncols = size(weight_matrix)
    rows = weight_matrix.rowval
    colptr = weight_matrix.colptr
    vals = weight_matrix.nzval
    essential_vals = zeros(T, length(vals))
    nonessential_vals = zeros(T, length(vals))

    num_inds = num_industries(info)
    essential_by_industry = zeros(T, num_inds)

    @inbounds for col = 1:ncols
        start_idx = colptr[col]
        stop_idx = colptr[col + 1] - 1
        if start_idx > stop_idx
            continue
        end

        fill!(essential_by_industry, zero(T))
        all_suppliers_total = zero(T)
        for idx = start_idx:stop_idx
            row = rows[idx]
            val = vals[idx]
            all_suppliers_total += val
            if is_essential(info, row)
                essential_by_industry[get_industry(info, row)] += val
            end
        end

        for idx = start_idx:stop_idx
            row = rows[idx]
            val = vals[idx]
            if is_essential(info, row)
                denom = essential_by_industry[get_industry(info, row)]
                essential_vals[idx] = denom == 0 ? zero(T) : val / denom
            else
                nonessential_vals[idx] = all_suppliers_total == 0 ? zero(T) : val / all_suppliers_total
            end
        end
    end

    essential_csc = SparseMatrixCSC(nrows, ncols, copy(colptr), copy(rows), essential_vals)
    nonessential_csc = SparseMatrixCSC(nrows, ncols, copy(colptr), copy(rows), nonessential_vals)

    I1, J1, V1 = findnz(essential_csc)
    I2, J2, V2 = findnz(nonessential_csc)
    return sparsecsr(I1, J1, V1, nrows, ncols), sparsecsr(I2, J2, V2, nrows, ncols)
end

function _validate_weight_matrix_entries(weight_matrix::AbstractMatrix{T}) where {T<:Real}
    @inbounds for i in eachindex(weight_matrix)
        v = weight_matrix[i]
        if !isfinite(v) || v < zero(T)
            throw(DomainError(v, "weight_matrix entries must be finite and nonnegative"))
        end
    end
    return nothing
end

function _validate_weight_matrix_entries(weight_matrix::SparseMatrixCSC{T}) where {T<:Real}
    @inbounds for i in eachindex(weight_matrix.nzval)
        v = weight_matrix.nzval[i]
        if !isfinite(v) || v < zero(T)
            throw(DomainError(v, "weight_matrix entries must be finite and nonnegative"))
        end
    end
    return nothing
end

"""
    ESRIEconomy(weight_matrix, info::IndustryInfo)

Precompute normalized upstream/downstream impact operators and firm output weights.
"""
function ESRIEconomy(weight_matrix::AbstractMatrix{T}, info::IndustryInfo) where {T<:Real}
    n = length(info)
    if size(weight_matrix, 1) != n || size(weight_matrix, 2) != n
        throw(DimensionMismatch("weight_matrix must be square and match info size"))
    end
    _validate_weight_matrix_entries(weight_matrix)

    upstream_impact = create_upstream_impact_matrix(weight_matrix)
    downstream_impact_essential, downstream_impact_nonessential = compute_downstream_impact_matrices(weight_matrix, info)

    column_sums = vec(sum(weight_matrix, dims = 1))
    row_sums = vec(sum(weight_matrix, dims = 2))
    total_output = sum(column_sums)

    return ESRIEconomy(
        info,
        upstream_impact,
        downstream_impact_essential,
        downstream_impact_nonessential,
        column_sums,
        row_sums,
        total_output,
        n,
    )
end
