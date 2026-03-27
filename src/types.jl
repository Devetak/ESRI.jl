struct IndustryInfo{TI<:AbstractVector{Int},TB<:AbstractVector{Bool}}
    industry_of_firm::TI
    essential_industry::TB
    essential_firm::TB
end

struct ESRIEconomy{T,TU,TD}
    info::IndustryInfo
    upstream_impact::TU
    downstream_impact::TD
    column_sums::Vector{T}
    row_sums::Vector{T}
    total_output::T
    n::Int
end

struct ESRIResult{T}
    esri::T
    upstream::Vector{T}
    downstream::Vector{T}
end

function IndustryInfo(
    industry_of_firm::AbstractVector{Int},
    essential_industry::AbstractVector{Bool},
)
    num_firms = length(industry_of_firm)
    num_industries = length(essential_industry)
    @assert num_industries > 0 "essential_industry must be non-empty"
    @assert all(1 .<= industry_of_firm .<= num_industries) "industry indices must be in 1:num_industries"
    essential_firm = similar(essential_industry, Bool, num_firms)
    @inbounds for i in eachindex(industry_of_firm)
        essential_firm[i] = essential_industry[industry_of_firm[i]]
    end
    return IndustryInfo(industry_of_firm, essential_industry, essential_firm)
end

function ESRIEconomy(weight_matrix::AbstractMatrix{T}, info::IndustryInfo) where {T<:Real}
    n = length(info)
    @assert size(weight_matrix, 1) == n && size(weight_matrix, 2) == n "weight_matrix must be square and match info size"

    upstream_impact = create_upstream_impact_matrix(weight_matrix)
    downstream_impact = if weight_matrix isa SparseMatrixCSC{T}
        compute_downstream_impact_matrix_csr(weight_matrix, info)
    else
        compute_downstream_impact_matrix(weight_matrix, info)
    end

    column_sums = vec(sum(weight_matrix, dims = 1))
    row_sums = vec(sum(weight_matrix, dims = 2))
    total_output = sum(column_sums)

    return ESRIEconomy(
        info,
        upstream_impact,
        downstream_impact,
        column_sums,
        row_sums,
        total_output,
        n,
    )
end

Base.length(info::IndustryInfo) = length(info.industry_of_firm)
Base.length(econ::ESRIEconomy) = econ.n

num_industries(info::IndustryInfo) = length(info.essential_industry)

@inline function is_essential(info::IndustryInfo, idx::Integer)
    @boundscheck checkbounds(info.essential_firm, idx)
    @inbounds return info.essential_firm[idx]
end

@inline function get_industry(info::IndustryInfo, idx::Integer)
    @boundscheck checkbounds(info.industry_of_firm, idx)
    @inbounds return info.industry_of_firm[idx]
end
