"""
    IndustryInfo(industry_of_firm, essential_industry)
    IndustryInfo(industry_of_firm, input_classification)
    IndustryInfo(industry_of_firm)

Firm industry ids and input classifications. `input_classification[supplier, customer]`
is `0` (no downstream impact), `1` (nonessential), or `2` (essential). The
legacy `essential_industry` and `essential_firm` fields remain for compatibility;
with matrix-based customer-specific classifications they are true only for rows
consisting entirely of `2`s. The one-argument form classifies every pair as `1`,
giving a purely linear baseline.
"""
struct IndustryInfo{TI<:AbstractVector{Int},TB<:AbstractVector{Bool}}
    industry_of_firm::TI
    essential_industry::TB
    essential_firm::TB
    input_classification::Matrix{UInt8}
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

struct _ESRIWorkspace{T}
    current_upstream::Vector{T}
    previous_upstream::Vector{T}
    current_downstream::Vector{T}
    previous_downstream::Vector{T}
    sigmas::Vector{T}
    essential_matrix::Matrix{T}
    essential_touched::Vector{Int}
    temp_sums::Vector{T}
    nonessential_vector::Vector{T}
    psi::Vector{T}
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

function IndustryInfo(industry_of_firm::AbstractVector{<:Integer})
    isempty(industry_of_firm) && throw(ArgumentError("industry_of_firm must be non-empty"))
    return IndustryInfo(industry_of_firm, falses(maximum(industry_of_firm)))
end

function _industry_indices(industry_of_firm::AbstractVector{<:Integer}, nindustries::Int)
    nindustries > 0 || throw(ArgumentError("essential_industry must be non-empty"))
    indices = Int.(industry_of_firm)
    all(idx -> 1 <= idx <= nindustries, indices) || throw(
        ArgumentError("industry_of_firm values must be in 1:length(essential_industry)"),
    )
    return indices
end

function IndustryInfo(
    industry_of_firm::AbstractVector{<:Integer},
    essential_industry::AbstractVector{Bool},
)
    nindustries = length(essential_industry)
    firm_industry = _industry_indices(industry_of_firm, nindustries)
    essential_industry = Vector{Bool}(essential_industry)
    return IndustryInfo(
        firm_industry,
        essential_industry,
        essential_industry[firm_industry],
    )
end

# Keep the pre-classification three-field constructor usable.  Its matrix view
# has the same row-constant 2/1 semantics as the old Boolean API.
function IndustryInfo{TI,TB}(
    industry_of_firm::TI,
    essential_industry::TB,
    essential_firm::TB,
) where {TI<:AbstractVector{Int},TB<:AbstractVector{Bool}}
    nindustries = length(essential_industry)
    return IndustryInfo{TI,TB}(
        industry_of_firm,
        essential_industry,
        essential_firm,
        repeat(UInt8.(essential_industry) .+ UInt8(1), 1, nindustries),
    )
end

IndustryInfo(
    industry_of_firm::TI,
    essential_industry::TB,
    essential_firm::TB,
) where {TI<:AbstractVector{Int},TB<:AbstractVector{Bool}} =
    IndustryInfo{TI,TB}(industry_of_firm, essential_industry, essential_firm)

function IndustryInfo(
    industry_of_firm::AbstractVector{<:Integer},
    input_classification::AbstractMatrix{Bool},
)
    throw(
        ArgumentError(
            "Boolean input_classification is ambiguous; use integer codes 0, 1, and 2",
        ),
    )
end

function IndustryInfo(
    industry_of_firm::AbstractVector{<:Integer},
    input_classification::AbstractMatrix{<:Integer},
)
    nindustries, ncustomers = size(input_classification)
    nindustries == ncustomers ||
        throw(DimensionMismatch("input_classification must be square"))
    classification = Matrix{UInt8}(input_classification)
    essential_industry = vec(all(==(UInt8(2)), classification; dims = 2))
    firm_industry = _industry_indices(industry_of_firm, nindustries)
    return IndustryInfo(
        firm_industry,
        essential_industry,
        essential_industry[firm_industry],
        classification,
    )
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
