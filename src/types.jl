"""
    IndustryInfo(industry_of_firm, essential_industry)

Firm industry ids and per-industry essentiality flags.
"""
struct IndustryInfo{TI<:AbstractVector{Int},TB<:AbstractVector{Bool}}
    industry_of_firm::TI
    essential_industry::TB
    essential_firm::TB
end

"""
    ESRIEconomy

Cached operators, weights, and totals for repeated ESRI runs.
"""
struct ESRIEconomy{T,I<:IndustryInfo,TU,TD,VT<:AbstractVector{T}}
    info::I
    upstream_impact::TU
    downstream_impact_essential::TD
    downstream_impact_nonessential::TD
    column_sums::VT
    row_sums::VT
    total_output::T
    n::Int
end

"""
    ESRIResult

Single-scenario ESRI plus converged upstream and downstream states.
"""
struct ESRIResult{T}
    esri::T
    upstream::Vector{T}
    downstream::Vector{T}
end

"""
    IndustryInfo(industry_of_firm::AbstractVector{<:Integer}, essential_industry::AbstractVector{Bool})

Build immutable industry metadata for ESRI.
`industry_of_firm[i]` is the 1-based industry id of firm `i`.
`essential_industry[k]` marks whether industry `k` is essential.
"""
function IndustryInfo(
    industry_of_firm::AbstractVector{<:Integer},
    essential_industry::AbstractVector{Bool},
)
    if isempty(essential_industry)
        throw(ArgumentError("essential_industry must be non-empty"))
    end

    firm_industry = Vector{Int}(undef, length(industry_of_firm))
    essential_industry_vec = Vector{Bool}(essential_industry)
    max_industry = length(essential_industry_vec)

    @inbounds for i in eachindex(industry_of_firm)
        idx = Int(industry_of_firm[i])
        if idx < 1 || idx > max_industry
            throw(ArgumentError("industry_of_firm values must be in 1:length(essential_industry)"))
        end
        firm_industry[i] = idx
    end

    essential_firm = Vector{Bool}(undef, length(firm_industry))
    @inbounds for i in eachindex(firm_industry)
        essential_firm[i] = essential_industry_vec[firm_industry[i]]
    end

    return IndustryInfo(firm_industry, essential_industry_vec, essential_firm)
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
