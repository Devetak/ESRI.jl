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

    # Build transposed sparse matrix: result[target, source] = weight[source, target] / row_sums[source]
    # We need to transpose, so we build a new CSC matrix with swapped dimensions
    nnz = length(vals)
    new_rows = Vector{Int}(undef, nnz)
    new_cols = Vector{Int}(undef, nnz)
    new_vals = Vector{T}(undef, nnz)

    idx_count = 0
    @inbounds for col = 1:ncols
        for idx = colptr[col]:(colptr[col+1]-1)
            source = rows[idx]  # source row in original
            target = col        # target col in original (becomes row in result)
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

    # Resize arrays if needed
    resize!(new_rows, idx_count)
    resize!(new_cols, idx_count)
    resize!(new_vals, idx_count)

    # Build CSC matrix from COO format (rows, cols, vals)
    # Result should be ncols x nrows (transposed)
    return sparse(new_rows, new_cols, new_vals, ncols, nrows)
end

function compute_downstream_impact_matrix(
    weight_matrix::AbstractMatrix{T},
    info::IndustryInfo,
) where {T<:Real}
    nrows, ncols = size(weight_matrix)
    result = zeros(T, nrows, ncols)
    column_sums = vec(sum(weight_matrix, dims = 1))
    num_inds = num_industries(info)
    partial = zeros(T, num_inds)

    @inbounds for col = 1:ncols
        col_sum = column_sums[col]
        if col_sum == 0
            continue
        end
        fill!(partial, zero(T))
        col_is_essential = is_essential(info, col)
        if col_is_essential
            for idx = 1:nrows
                val = weight_matrix[idx, col]
                if val == 0
                    continue
                end
                if is_essential(info, idx)
                    partial[get_industry(info, idx)] += val
                end
            end
        end

        for idx = 1:nrows
            val = weight_matrix[idx, col]
            if val == 0
                continue
            end
            if col_is_essential && is_essential(info, idx)
                denom = partial[get_industry(info, idx)]
                result[idx, col] = denom == 0 ? zero(T) : val / denom
            else
                result[idx, col] = val / col_sum
            end
        end
    end
    return result
end

function compute_downstream_impact_matrix(
    weight_matrix::SparseMatrixCSC{T},
    info::IndustryInfo,
) where {T<:Real}
    nrows, ncols = size(weight_matrix)
    rows = weight_matrix.rowval
    colptr = weight_matrix.colptr
    vals = weight_matrix.nzval
    result_vals = similar(vals)

    num_inds = num_industries(info)
    partial = zeros(T, num_inds)

    @inbounds for col = 1:ncols
        start_idx = colptr[col]
        stop_idx = colptr[col+1] - 1
        if start_idx > stop_idx
            continue
        end
        col_sum = zero(T)
        for idx = start_idx:stop_idx
            col_sum += vals[idx]
        end
        if col_sum == 0
            for idx = start_idx:stop_idx
                result_vals[idx] = zero(T)
            end
            continue
        end

        col_is_essential = is_essential(info, col)
        if col_is_essential
            fill!(partial, zero(T))
            for idx = start_idx:stop_idx
                row = rows[idx]
                if is_essential(info, row)
                    partial[get_industry(info, row)] += vals[idx]
                end
            end
        end

        for idx = start_idx:stop_idx
            row = rows[idx]
            val = vals[idx]
            if col_is_essential && is_essential(info, row)
                denom = partial[get_industry(info, row)]
                result_vals[idx] = denom == 0 ? zero(T) : val / denom
            else
                result_vals[idx] = val / col_sum
            end
        end
    end

    return SparseMatrixCSC(nrows, ncols, copy(colptr), copy(rows), result_vals)
end

function compute_downstream_impact_matrix_csr(
    weight_matrix::SparseMatrixCSC{T},
    info::IndustryInfo,
) where {T<:Real}
    impact_csc = compute_downstream_impact_matrix(weight_matrix, info)
    # Convert CSC to COO format for sparsecsr
    I, J, V = findnz(impact_csc)
    return sparsecsr(I, J, V, size(impact_csc)...)
end
