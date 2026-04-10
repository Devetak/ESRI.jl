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
        if row_sums[firm_idx] == 0
            sigmas[firm_idx] = zeroT
        elseif denom == 0
            sigmas[firm_idx] = oneT
        else
            sigmas[firm_idx] = min(row_sums[firm_idx] / denom, oneT)
        end
    end
    return sigmas
end

function _accumulate_downstream_components!(
    essential_matrix::AbstractMatrix,
    nonessential_vector::AbstractVector,
    downstream_vector::AbstractVector,
    sigmas::AbstractVector,
    essential_impact::AbstractMatrix,
    nonessential_impact::AbstractMatrix,
    info::IndustryInfo,
    industry_map::AbstractVector{Int} = info.industry_of_firm,
)
    T = eltype(essential_matrix)
    oneT = one(T)

    fill!(essential_matrix, zero(T))
    fill!(nonessential_vector, zero(T))

    nfirms = length(downstream_vector)
    ncustomers = size(essential_matrix, 1)
    @inbounds for supplier = 1:nfirms
        industry = industry_map[supplier]
        factor = sigmas[supplier] * (oneT - downstream_vector[supplier])
        if factor == 0
            continue
        end
        @simd for customer = 1:ncustomers
            essential_matrix[customer, industry] += factor * essential_impact[supplier, customer]
            nonessential_vector[customer] += factor * nonessential_impact[supplier, customer]
        end
    end
    return nothing
end

function _accumulate_downstream_components!(
    essential_matrix::AbstractMatrix,
    nonessential_vector::AbstractVector,
    downstream_vector::AbstractVector,
    sigmas::AbstractVector,
    essential_impact::SparseMatrixCSR,
    nonessential_impact::SparseMatrixCSR,
    info::IndustryInfo,
    industry_map::AbstractVector{Int} = info.industry_of_firm,
)
    T = eltype(essential_matrix)
    oneT = one(T)

    fill!(essential_matrix, zero(T))
    fill!(nonessential_vector, zero(T))

    essential_rowptr = essential_impact.rowptr
    essential_colval = essential_impact.colval
    essential_vals = essential_impact.nzval

    nonessential_rowptr = nonessential_impact.rowptr
    nonessential_colval = nonessential_impact.colval
    nonessential_vals = nonessential_impact.nzval

    nfirms = length(downstream_vector)
    @inbounds for supplier = 1:nfirms
        industry = industry_map[supplier]
        factor = sigmas[supplier] * (oneT - downstream_vector[supplier])
        if factor == 0
            continue
        end
        for idx = essential_rowptr[supplier]:(essential_rowptr[supplier + 1] - 1)
            customer = essential_colval[idx]
            essential_matrix[customer, industry] += factor * essential_vals[idx]
        end
        for idx = nonessential_rowptr[supplier]:(nonessential_rowptr[supplier + 1] - 1)
            customer = nonessential_colval[idx]
            nonessential_vector[customer] += factor * nonessential_vals[idx]
        end
    end
    return nothing
end

function downstream_step!(
    downstream_vector::AbstractVector,
    essential_matrix::AbstractMatrix,
    nonessential_vector::AbstractVector,
    psi::AbstractVector,
)
    n = size(essential_matrix, 1)
    ncols = size(essential_matrix, 2)
    @inbounds for firm_idx = 1:n
        essential_shortage = zero(eltype(downstream_vector))
        for col = 1:ncols
            val = essential_matrix[firm_idx, col]
            if val > essential_shortage
                essential_shortage = val
            end
        end
        essential_health = one(eltype(downstream_vector)) - essential_shortage
        nonessential_health = one(eltype(downstream_vector)) - nonessential_vector[firm_idx]
        downstream_vector[firm_idx] = min(essential_health, nonessential_health, psi[firm_idx])
    end
    return downstream_vector
end

function downstream_shock!(
    essential_impact,
    nonessential_impact,
    info::IndustryInfo,
    row_sums::AbstractVector,
    psi::AbstractVector,
    sigmas::AbstractVector,
    essential_matrix::AbstractMatrix,
    nonessential_vector::AbstractVector,
    current_downstream::AbstractVector,
    temp_sums::AbstractVector;
    industry_map::AbstractVector{Int} = info.industry_of_firm,
)
    compute_sigmas!(
        sigmas,
        row_sums,
        current_downstream,
        info,
        temp_sums,
        industry_map,
    )

    _accumulate_downstream_components!(
        essential_matrix,
        nonessential_vector,
        current_downstream,
        sigmas,
        essential_impact,
        nonessential_impact,
        info,
        industry_map,
    )
    downstream_step!(current_downstream, essential_matrix, nonessential_vector, psi)
    return nothing
end
