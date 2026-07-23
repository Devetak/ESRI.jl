"""
    IndustryInfo(industry_of_firm, essential_industry)
    IndustryInfo(industry_of_firm, input_classification)
    IndustryInfo(industry_of_firm)

Firm industry ids and input classifications. `input_classification[supplier, customer]`
is `0` (no downstream impact), `1` (nonessential), or `2` (essential).
"""
struct IndustryInfo{TI<:AbstractVector{Int},TB<:AbstractVector{Bool},TC<:AbstractMatrix{UInt8}}
    industry_of_firm::TI
    essential_industry::TB
    essential_firm::TB
    input_classification::TC
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
    IndustryInfo(industry_of_firm::AbstractVector{<:Integer})

Build immutable industry metadata with every supplier-customer industry pair
classified as essential. The number of industries is inferred from the largest
industry id.
"""
function IndustryInfo(industry_of_firm::AbstractVector{<:Integer})
    isempty(industry_of_firm) && throw(ArgumentError("industry_of_firm must be non-empty"))
    nindustries = maximum(industry_of_firm)
    nindustries > 0 || throw(ArgumentError("industry_of_firm values must be positive"))
    return IndustryInfo(industry_of_firm, trues(nindustries))
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

    nindustries = length(essential_industry)
    input_classification = Matrix{UInt8}(undef, nindustries, nindustries)
    @inbounds for supplier_industry in 1:nindustries
        fill!(
            @view(input_classification[supplier_industry, :]),
            essential_industry[supplier_industry] ? UInt8(2) : UInt8(1),
        )
    end
    return IndustryInfo(industry_of_firm, input_classification)
end

"""
    IndustryInfo(industry_of_firm::AbstractVector{<:Integer}, input_classification::AbstractMatrix{<:Integer})

Build immutable industry metadata with a supplier-industry by customer-industry
classification matrix. Entries must be `0`, `1`, or `2`.
"""
function IndustryInfo(
    industry_of_firm::AbstractVector{<:Integer},
    input_classification::AbstractMatrix{<:Integer},
)
    nindustries, ncustomers = size(input_classification)
    nindustries == ncustomers || throw(DimensionMismatch("input_classification must be square"))
    nindustries > 0 || throw(ArgumentError("input_classification must be non-empty"))

    firm_industry = Vector{Int}(undef, length(industry_of_firm))
    @inbounds for i in eachindex(industry_of_firm)
        idx = Int(industry_of_firm[i])
        if idx < 1 || idx > nindustries
            throw(ArgumentError("industry_of_firm values must be in 1:size(input_classification, 1)"))
        end
        firm_industry[i] = idx
    end

    classification = Matrix{UInt8}(undef, nindustries, nindustries)
    essential_industry_vec = Vector{Bool}(undef, nindustries)
    @inbounds for supplier_industry in 1:nindustries
        universally_essential = true
        for customer_industry in 1:nindustries
            value = input_classification[supplier_industry, customer_industry]
            if value != 0 && value != 1 && value != 2
                throw(ArgumentError("input_classification entries must be 0, 1, or 2"))
            end
            classification[supplier_industry, customer_industry] = UInt8(value)
            universally_essential &= value == 2
        end
        essential_industry_vec[supplier_industry] = universally_essential
    end

    essential_firm = Vector{Bool}(undef, length(firm_industry))
    @inbounds for i in eachindex(firm_industry)
        essential_firm[i] = essential_industry_vec[firm_industry[i]]
    end

    return IndustryInfo(firm_industry, essential_industry_vec, essential_firm, classification)
end

Base.length(info::IndustryInfo) = length(info.industry_of_firm)
Base.length(econ::ESRIEconomy) = econ.n

num_industries(info::IndustryInfo) = size(info.input_classification, 1)

@inline function is_essential(info::IndustryInfo, idx::Integer)
    @boundscheck checkbounds(info.essential_firm, idx)
    @inbounds return info.essential_firm[idx]
end

@inline function get_industry(info::IndustryInfo, idx::Integer)
    @boundscheck checkbounds(info.industry_of_firm, idx)
    @inbounds return info.industry_of_firm[idx]
end
