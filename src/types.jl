struct IndustryInfo{TI<:AbstractVector{Int},TB<:AbstractVector{Bool}}
    industry_of_firm::TI
    essential_industry::TB
    essential_firm::TB
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

Base.length(info::IndustryInfo) = length(info.industry_of_firm)

num_industries(info::IndustryInfo) = length(info.essential_industry)

@inline function is_essential(info::IndustryInfo, idx::Integer)
    @boundscheck checkbounds(info.essential_firm, idx)
    @inbounds return info.essential_firm[idx]
end

@inline function get_industry(info::IndustryInfo, idx::Integer)
    @boundscheck checkbounds(info.industry_of_firm, idx)
    @inbounds return info.industry_of_firm[idx]
end
