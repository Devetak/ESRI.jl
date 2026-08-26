const _IHS_DATA_DIR = joinpath(@__DIR__, "..", "data")
const _IHS_CLASSIFICATION_CACHE = Ref{Union{Nothing,Matrix{UInt8}}}(nothing)

"""
    ihs_input_classification()

Return a copy of the bundled 616 × 616 IHS supplier-by-customer input-classification
matrix. Entries are `0`, `1`, and `2`. Use [`ihs_industry_codes`](@ref) for the
row and column order.
"""
function ihs_input_classification()
    cached = _IHS_CLASSIFICATION_CACHE[]
    cached !== nothing && return copy(cached)

    bytes = read(joinpath(_IHS_DATA_DIR, "ihs_classification.bin"))
    n = Int(bytes[1]) + (Int(bytes[2]) << 8)
    profiles = Int(bytes[3]) + (Int(bytes[4]) << 8)
    length(bytes) == 4 + n + n * profiles || error("invalid bundled IHS classification")

    classification = Matrix{UInt8}(undef, n, n)
    @inbounds for customer = 1:n
        profile = Int(bytes[4+customer])
        source = 5 + n + profile * n
        copyto!(classification, (customer - 1) * n + 1, bytes, source, n)
    end
    _IHS_CLASSIFICATION_CACHE[] = classification
    return copy(classification)
end

"""
    ihs_industry_codes()

Return the industry-code order used by [`ihs_input_classification`](@ref).
"""
ihs_industry_codes() = readlines(joinpath(_IHS_DATA_DIR, "ihs_industry_codes.txt"))
