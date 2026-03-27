function compute_sigmas!(
    sigmas::AbstractVector,
    row_sums::AbstractVector,
    downstream_vector::AbstractVector,
    info::IndustryInfo,
    temp_sums::AbstractVector,
    industry_map::AbstractVector{Int} = info.industry_of_firm,
)
    T = eltype(sigmas)
    zeroT = zero(T)
    oneT = one(T)

    fill!(temp_sums, zeroT)
    @inbounds for firm_idx in eachindex(sigmas)
        industry = industry_map[firm_idx]
        temp_sums[industry] += row_sums[firm_idx] * downstream_vector[firm_idx]
    end
    @inbounds for firm_idx in eachindex(sigmas)
        industry = industry_map[firm_idx]
        denom = temp_sums[industry]
        if denom == 0
            sigmas[firm_idx] = oneT
        else
            sigmas[firm_idx] = min(row_sums[firm_idx] / denom, oneT)
        end
    end
    return sigmas
end

function compute_product_matrix!(
    product_matrix::AbstractMatrix,
    downstream_vector::AbstractVector,
    sigmas::AbstractVector,
    impact_matrix::AbstractMatrix,
    info::IndustryInfo,
    industry_map::AbstractVector{Int} = info.industry_of_firm,
    one_minus_downstream::AbstractVector = similar(downstream_vector),
)
    T = eltype(product_matrix)
    oneT = one(T)

    fill!(product_matrix, oneT)

    @assert length(one_minus_downstream) == length(downstream_vector) "one_minus_downstream length must match downstream_vector"
    @inbounds for i in eachindex(downstream_vector)
        one_minus_downstream[i] = oneT - downstream_vector[i]
    end

    nfirms = length(downstream_vector)
    @inbounds for col = 1:nfirms
        industry = industry_map[col]
        factor = sigmas[col] * one_minus_downstream[col]
        if factor == 0
            continue
        end
        for row = 1:size(product_matrix, 1)
            val = impact_matrix[col, row]
            if val != 0
                product_matrix[row, industry] -= factor * val
            end
        end
    end
    return product_matrix
end

function compute_product_matrix!(
    product_matrix::AbstractMatrix,
    downstream_vector::AbstractVector,
    sigmas::AbstractVector,
    impact_matrix::SparseMatrixCSC,
    info::IndustryInfo,
    industry_map::AbstractVector{Int} = info.industry_of_firm,
    one_minus_downstream::AbstractVector = similar(downstream_vector),
)
    T = eltype(product_matrix)
    oneT = one(T)

    fill!(product_matrix, oneT)

    @assert length(one_minus_downstream) == length(downstream_vector) "one_minus_downstream length must match downstream_vector"
    @inbounds for i in eachindex(downstream_vector)
        one_minus_downstream[i] = oneT - downstream_vector[i]
    end

    rows = impact_matrix.rowval
    colptr = impact_matrix.colptr
    vals = impact_matrix.nzval

    ncols = size(impact_matrix, 2)
    @inbounds for col = 1:ncols
        industry = industry_map[col]
        factor = sigmas[col] * one_minus_downstream[col]
        if factor == 0
            continue
        end
        for idx = colptr[col]:(colptr[col+1]-1)
            row = rows[idx]
            product_matrix[row, industry] -= factor * vals[idx]
        end
    end
    return product_matrix
end

function compute_product_matrix!(
    product_matrix::AbstractMatrix,
    downstream_vector::AbstractVector,
    sigmas::AbstractVector,
    impact_matrix::SparseMatrixCSR,
    info::IndustryInfo,
    industry_map::AbstractVector{Int} = info.industry_of_firm,
    one_minus_downstream::AbstractVector = similar(downstream_vector),
)
    T = eltype(product_matrix)
    oneT = one(T)

    fill!(product_matrix, oneT)

    @assert length(one_minus_downstream) == length(downstream_vector) "one_minus_downstream length must match downstream_vector"
    @inbounds for i in eachindex(downstream_vector)
        one_minus_downstream[i] = oneT - downstream_vector[i]
    end

    rowptr = impact_matrix.rowptr
    colval = impact_matrix.colval
    vals = impact_matrix.nzval

    nrows = size(impact_matrix, 1)
    @inbounds for row = 1:nrows
        industry = industry_map[row]
        factor = sigmas[row] * one_minus_downstream[row]
        if factor == 0
            continue
        end
        for idx = rowptr[row]:(rowptr[row+1]-1)
            col = colval[idx]
            product_matrix[col, industry] -= factor * vals[idx]
        end
    end
    return product_matrix
end

function downstream_step!(
    downstream_vector::AbstractVector,
    product_matrix::AbstractMatrix,
    info::IndustryInfo,
    essential_mask::AbstractVector{Bool} = info.essential_firm,
)
    n = size(product_matrix, 1)
    ncols = size(product_matrix, 2)
    @inbounds for firm_idx = 1:n
        if essential_mask[firm_idx]
            # Compute minimum manually to avoid view allocation
            min_val = product_matrix[firm_idx, 1]
            for col = 2:ncols
                val = product_matrix[firm_idx, col]
                if val < min_val
                    min_val = val
                end
            end
            downstream_vector[firm_idx] = min_val
        else
            downstream_vector[firm_idx] = product_matrix[firm_idx, end]
        end
    end
    return downstream_vector
end

function downstream_shock!(
    impact_matrix,
    info::IndustryInfo,
    row_sums::AbstractVector,
    firm_idx::Integer,
    sigmas::AbstractVector,
    product_matrix::AbstractMatrix,
    current_downstream::AbstractVector,
    previous_downstream::AbstractVector,
    temp_sums::AbstractVector;
    maxiter::Int = 100,
    tol::Real = 1e-2,
    verbose::Bool = false,
    every::Int = 10,
)
    T = eltype(current_downstream)
    zeroT = zero(T)
    oneT = one(T)

    # Precompute industry_map and essential_mask
    industry_map = info.industry_of_firm
    essential_mask = info.essential_firm

    # Preallocate one_minus_downstream buffer
    one_minus_downstream = similar(current_downstream)

    current_downstream .= oneT
    previous_downstream .= oneT
    current_downstream[firm_idx] = zeroT
    previous_downstream[firm_idx] = zeroT

    for iter = 1:maxiter
        compute_sigmas!(
            sigmas,
            row_sums,
            current_downstream,
            info,
            temp_sums,
            industry_map,
        )

        compute_product_matrix!(
            product_matrix,
            current_downstream,
            sigmas,
            impact_matrix,
            info,
            industry_map,
            one_minus_downstream,
        )
        downstream_step!(current_downstream, product_matrix, info, essential_mask)
        current_downstream[firm_idx] = zeroT

        distance = _linf_distance(current_downstream, previous_downstream)
        if distance < tol
            return nothing
        end
        if verbose && ((iter - 1) % every == 0)
            @info "downstream iteration" iter distance
        end
        copyto!(previous_downstream, current_downstream)
    end
    return nothing
end
