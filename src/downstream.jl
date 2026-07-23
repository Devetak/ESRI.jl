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

function _accumulate_active_sparse_downstream_components!(
    essential_matrix::Matrix{T},
    essential_touched::Vector{Int},
    nonessential_vector::AbstractVector{T},
    downstream_vector::AbstractVector{T},
    sigmas::AbstractVector{T},
    essential_impact::SparseMatrixCSR,
    nonessential_impact::SparseMatrixCSR,
    info::IndustryInfo,
    industry_map::AbstractVector{Int} = info.industry_of_firm,
) where {T}
    zeroT = zero(T)
    oneT = one(T)

    fill!(nonessential_vector, zeroT)

    essential_rowptr = essential_impact.rowptr
    essential_colval = essential_impact.colval
    essential_vals = essential_impact.nzval
    nonessential_rowptr = nonessential_impact.rowptr
    nonessential_colval = nonessential_impact.colval
    nonessential_vals = nonessential_impact.nzval

    sizehint!(essential_touched, length(essential_vals))
    nfirms = length(downstream_vector)
    ncustomers = size(essential_matrix, 1)
    @inbounds for supplier = 1:nfirms
        factor = sigmas[supplier] * (oneT - downstream_vector[supplier])
        factor == zeroT && continue
        industry = industry_map[supplier]

        for idx = essential_rowptr[supplier]:(essential_rowptr[supplier + 1] - 1)
            linear_idx = essential_colval[idx] + (industry - 1) * ncustomers
            if essential_matrix[linear_idx] == zeroT
                push!(essential_touched, linear_idx)
            end
            essential_matrix[linear_idx] += factor * essential_vals[idx]
        end

        for idx = nonessential_rowptr[supplier]:(nonessential_rowptr[supplier + 1] - 1)
            nonessential_vector[nonessential_colval[idx]] += factor * nonessential_vals[idx]
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
    fill!(downstream_vector, zero(eltype(downstream_vector)))
    @inbounds for col = 1:ncols
        @simd for firm_idx = 1:n
            downstream_vector[firm_idx] = max(
                downstream_vector[firm_idx],
                essential_matrix[firm_idx, col],
            )
        end
    end
    @inbounds for firm_idx = 1:n
        essential_health = one(eltype(downstream_vector)) - downstream_vector[firm_idx]
        nonessential_health = one(eltype(downstream_vector)) - nonessential_vector[firm_idx]
        downstream_vector[firm_idx] = min(essential_health, nonessential_health, psi[firm_idx])
    end
    return downstream_vector
end

function downstream_step!(
    downstream_vector::AbstractVector{T},
    essential_matrix::Matrix{T},
    essential_touched::Vector{Int},
    nonessential_vector::AbstractVector{T},
    psi::AbstractVector{T},
) where {T}
    zeroT = zero(T)
    oneT = one(T)
    n = size(essential_matrix, 1)
    fill!(downstream_vector, zeroT)
    @inbounds for linear_idx in essential_touched
        customer = mod(linear_idx - 1, n) + 1
        downstream_vector[customer] = max(downstream_vector[customer], essential_matrix[linear_idx])
        essential_matrix[linear_idx] = zeroT
    end
    empty!(essential_touched)

    @inbounds for firm_idx in eachindex(downstream_vector)
        essential_health = oneT - downstream_vector[firm_idx]
        nonessential_health = oneT - nonessential_vector[firm_idx]
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

function downstream_shock!(
    essential_impact::SparseMatrixCSR,
    nonessential_impact::SparseMatrixCSR,
    info::IndustryInfo,
    row_sums::AbstractVector{T},
    psi::AbstractVector{T},
    sigmas::AbstractVector{T},
    essential_matrix::Matrix{T},
    nonessential_vector::AbstractVector{T},
    current_downstream::AbstractVector{T},
    temp_sums::AbstractVector{T};
    industry_map::AbstractVector{Int} = info.industry_of_firm,
) where {T}
    fill!(essential_matrix, zero(T))
    return downstream_shock!(
        essential_impact,
        nonessential_impact,
        info,
        row_sums,
        psi,
        sigmas,
        essential_matrix,
        nonessential_vector,
        current_downstream,
        temp_sums,
        Int[];
        industry_map = industry_map,
    )
end

function downstream_shock!(
    essential_impact,
    nonessential_impact,
    info::IndustryInfo,
    row_sums::AbstractVector,
    psi::AbstractVector,
    sigmas::AbstractVector,
    essential_matrix::Matrix,
    nonessential_vector::AbstractVector,
    current_downstream::AbstractVector,
    temp_sums::AbstractVector,
    essential_touched::Vector{Int};
    industry_map::AbstractVector{Int} = info.industry_of_firm,
)
    return downstream_shock!(
        essential_impact,
        nonessential_impact,
        info,
        row_sums,
        psi,
        sigmas,
        essential_matrix,
        nonessential_vector,
        current_downstream,
        temp_sums;
        industry_map = industry_map,
    )
end

function downstream_shock!(
    essential_impact::SparseMatrixCSR,
    nonessential_impact::SparseMatrixCSR,
    info::IndustryInfo,
    row_sums::AbstractVector{T},
    psi::AbstractVector{T},
    sigmas::AbstractVector{T},
    essential_matrix::Matrix{T},
    nonessential_vector::AbstractVector{T},
    current_downstream::AbstractVector{T},
    temp_sums::AbstractVector{T},
    essential_touched::Vector{Int};
    industry_map::AbstractVector{Int} = info.industry_of_firm,
) where {T}
    compute_sigmas!(sigmas, row_sums, current_downstream, info, temp_sums, industry_map)
    _accumulate_active_sparse_downstream_components!(
        essential_matrix,
        essential_touched,
        nonessential_vector,
        current_downstream,
        sigmas,
        essential_impact,
        nonessential_impact,
        info,
        industry_map,
    )
    downstream_step!(
        current_downstream,
        essential_matrix,
        essential_touched,
        nonessential_vector,
        psi,
    )
    return nothing
end

function _downstream_shock!(econ::ESRIEconomy{T}, workspace::_ESRIWorkspace{T}) where {T}
    return downstream_shock!(
        econ.downstream_impact_essential,
        econ.downstream_impact_nonessential,
        econ.info,
        econ.row_sums,
        workspace.psi,
        workspace.sigmas,
        workspace.essential_matrix,
        workspace.nonessential_vector,
        workspace.current_downstream,
        workspace.temp_sums,
        workspace.essential_touched,
    )
end
