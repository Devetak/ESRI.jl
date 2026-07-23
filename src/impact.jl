"""
    create_upstream_impact_matrix(weight_matrix)

Build the normalized upstream impact matrix.
"""
function create_upstream_impact_matrix(weight_matrix::AbstractMatrix{T}) where {T<:Real}
    nrows, ncols = size(weight_matrix)
    TF = float(T)
    result = zeros(TF, nrows, ncols)

    @inbounds for source = 1:nrows
        denom = zero(TF)
        for target = 1:ncols
            denom += weight_matrix[source, target]
        end
        denom == 0 && continue
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
    TF = float(T)
    nrows, ncols = size(weight_matrix)
    rows = weight_matrix.rowval
    colptr = weight_matrix.colptr
    vals = weight_matrix.nzval

    row_sums = zeros(TF, nrows)
    @inbounds for col = 1:ncols
        for idx = colptr[col]:(colptr[col+1]-1)
            row = rows[idx]
            row_sums[row] += vals[idx]
        end
    end

    nnz = length(vals)
    new_rows = Vector{Int}(undef, nnz)
    new_cols = Vector{Int}(undef, nnz)
    new_vals = Vector{TF}(undef, nnz)

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
    TF = float(T)
    nrows, ncols = size(weight_matrix)
    essential = zeros(TF, nrows, ncols)
    nonessential = zeros(TF, nrows, ncols)
    industry_of_firm = info.industry_of_firm
    input_classification = info.input_classification
    essential_by_industry = zeros(TF, num_industries(info))
    zeroT = zero(TF)

    @inbounds for col = 1:ncols
        fill!(essential_by_industry, zeroT)
        total = zeroT
        customer_industry = industry_of_firm[col]
        for row = 1:nrows
            val = weight_matrix[row, col]
            if val == 0
                continue
            end
            total += val
            supplier_industry = industry_of_firm[row]
            if input_classification[supplier_industry, customer_industry] == 2
                essential_by_industry[supplier_industry] += val
            end
        end

        for row = 1:nrows
            val = weight_matrix[row, col]
            if val == 0
                continue
            end
            supplier_industry = industry_of_firm[row]
            classification = input_classification[supplier_industry, customer_industry]
            if classification == 2
                denom = essential_by_industry[supplier_industry]
                essential[row, col] = denom == 0 ? zeroT : val / denom
            elseif classification == 1
                nonessential[row, col] = total == 0 ? zeroT : val / total
            end
        end
    end

    return essential, nonessential
end

function compute_downstream_impact_matrices(
    weight_matrix::SparseMatrixCSC{T},
    info::IndustryInfo,
) where {T<:Real}
    TF = float(T)
    nrows, ncols = size(weight_matrix)
    rows = weight_matrix.rowval
    colptr = weight_matrix.colptr
    vals = weight_matrix.nzval
    essential_vals = zeros(TF, length(vals))
    nonessential_vals = zeros(TF, length(vals))

    num_inds = num_industries(info)
    industry_of_firm = info.industry_of_firm
    input_classification = info.input_classification
    essential_by_industry = zeros(TF, num_inds)
    zeroT = zero(TF)

    @inbounds for col = 1:ncols
        start_idx = colptr[col]
        stop_idx = colptr[col + 1] - 1
        if start_idx > stop_idx
            continue
        end

        fill!(essential_by_industry, zeroT)
        all_suppliers_total = zeroT
        customer_industry = industry_of_firm[col]
        for idx = start_idx:stop_idx
            row = rows[idx]
            val = vals[idx]
            all_suppliers_total += val
            supplier_industry = industry_of_firm[row]
            if input_classification[supplier_industry, customer_industry] == 2
                essential_by_industry[supplier_industry] += val
            end
        end

        for idx = start_idx:stop_idx
            row = rows[idx]
            val = vals[idx]
            supplier_industry = industry_of_firm[row]
            classification = input_classification[supplier_industry, customer_industry]
            if classification == 2
                denom = essential_by_industry[supplier_industry]
                essential_vals[idx] = denom == 0 ? zeroT : val / denom
            elseif classification == 1
                nonessential_vals[idx] = all_suppliers_total == 0 ? zero(TF) : val / all_suppliers_total
            end
        end
    end

    essential_csc = SparseMatrixCSC(nrows, ncols, copy(colptr), copy(rows), essential_vals)
    nonessential_csc = SparseMatrixCSC(nrows, ncols, copy(colptr), copy(rows), nonessential_vals)
    dropzeros!(essential_csc)
    dropzeros!(nonessential_csc)
    return sparsecsr(essential_csc), sparsecsr(nonessential_csc)
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

function _promote_weight_matrix(weight_matrix::SparseMatrixCSC{T,Ti}) where {T<:Real,Ti<:Integer}
    TF = float(T)
    if T === TF
        return weight_matrix
    end
    return SparseMatrixCSC(
        size(weight_matrix, 1),
        size(weight_matrix, 2),
        copy(weight_matrix.colptr),
        copy(weight_matrix.rowval),
        TF.(weight_matrix.nzval),
    )
end

function _promote_weight_matrix(weight_matrix::AbstractMatrix{T}) where {T<:Real}
    TF = float(T)
    if T === TF
        return weight_matrix
    end
    return TF.(weight_matrix)
end

"""
    ESRIEconomy(weight_matrix, info::IndustryInfo)

Cache normalized upstream/downstream operators, output weights, and totals.
"""
function ESRIEconomy(weight_matrix::AbstractMatrix{T}, info::IndustryInfo) where {T<:Real}
    n = length(info)
    if size(weight_matrix, 1) != n || size(weight_matrix, 2) != n
        throw(DimensionMismatch("weight_matrix must be square and match info size"))
    end
    matrix = _promote_weight_matrix(weight_matrix)
    _validate_weight_matrix_entries(matrix)

    upstream_impact = create_upstream_impact_matrix(matrix)
    downstream_impact_essential, downstream_impact_nonessential = compute_downstream_impact_matrices(matrix, info)

    column_sums = vec(sum(matrix, dims = 1))
    row_sums = vec(sum(matrix, dims = 2))
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
