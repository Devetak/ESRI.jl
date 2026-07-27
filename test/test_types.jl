@testset "IndustryInfo" begin
    industry_ids = [1, 2, 1, 2]
    essential_industry = [true, false]
    info = IndustryInfo(industry_ids, essential_industry)

    @test ESRIcascade.num_industries(info) == 2
    @test info.industry_of_firm == industry_ids
    @test info.essential_industry == essential_industry
    @test info.essential_firm == [true, false, true, false]
    @test info isa IndustryInfo{Vector{Int},Vector{Bool}}

    legacy_info = IndustryInfo(industry_ids, essential_industry, [true, false, true, false])
    @test legacy_info isa IndustryInfo{Vector{Int},Vector{Bool}}
    @test legacy_info.input_classification == UInt8[2 2; 1 1]
    @test IndustryInfo{Vector{Int},Vector{Bool}}(
        industry_ids,
        essential_industry,
        [true, false, true, false],
    ).input_classification == legacy_info.input_classification

    default_info = IndustryInfo(industry_ids)
    @test default_info.input_classification == fill(UInt8(2), 2, 2)
    @test default_info.essential_firm == trues(length(industry_ids))
    @test_throws ArgumentError IndustryInfo(Int[])

    @test IndustryInfo(Int32[1, 2, 1, 2], Bool[true, false]).industry_of_firm ==
          [1, 2, 1, 2]
    @test_throws ArgumentError IndustryInfo([0, 1], [true, false])
    @test_throws ArgumentError IndustryInfo([3, 1], [true, false])
    @test_throws ArgumentError IndustryInfo([1, 1], Bool[])

    classification = [2 1; 0 2]
    classified = IndustryInfo(industry_ids, classification)
    @test classified.industry_of_firm == industry_ids
    @test classified.input_classification == UInt8[2 1; 0 2]
    @test ESRIcascade.num_industries(classified) == 2

    @test_throws DimensionMismatch IndustryInfo([1, 2], zeros(Int, 2, 3))
    @test_throws ArgumentError IndustryInfo([1, 2], zeros(Int, 0, 0))
    @test_throws ArgumentError IndustryInfo([3, 1], zeros(Int, 2, 2))
    @test_throws ArgumentError IndustryInfo([1, 2], Bool[true false; false true])
end

@testset "Bundled IHS classification" begin
    classification = ihs_input_classification()
    data_dir = joinpath(@__DIR__, "..", "data")
    @test size(classification) == (616, 616)
    @test ihs_industry_codes()[[1, end]] == ["0111", "9999"]
    @test count(==(UInt8(0)), classification) == 269204
    @test count(==(UInt8(1)), classification) == 45013
    @test count(==(UInt8(2)), classification) == 65239
    classification[1, 1] = 0
    @test ihs_input_classification()[1, 1] == 2
end

@testset "Bundled IHS C++ reference" begin
    classification = ihs_input_classification()
    industries = size(classification, 1)
    n = 2 * industries
    industry_of_firm = repeat(1:industries, inner = 2)
    closure_industries = round.(Int, range(1, industries; length = 16))
    closure_firms = 2 .* closure_industries .- 1
    @test industry_of_firm[closure_firms] == closure_industries

    large_n = 10_000
    firms_per_industry, extra_firms = divrem(large_n, industries)
    large_industry_of_firm =
        vcat(repeat(1:industries, inner = firms_per_industry), 1:extra_firms)
    large_closure_firms =
        [(industry - 1) * firms_per_industry + 1 for industry in closure_industries]
    @test large_industry_of_firm[large_closure_firms] == closure_industries
    @test length(unique(large_industry_of_firm[large_closure_firms])) == 16

    offsets = (1, 37, 91, 173, 311, 509, 701, 997)
    supplier = repeat(1:n, inner = length(offsets))
    customer = [mod1(i + offset, n) for i = 1:n for offset in offsets]
    weights = Float64.(1 .+ mod.(13 .* supplier .+ 7 .* customer, 17))
    econ = ESRIEconomy(
        sparse(supplier, customer, weights, n, n),
        IndustryInfo(industry_of_firm, classification),
    )
    scores = esri(
        econ;
        firm_indices = closure_firms,
        maxiter = 300,
        tol = 1e-2,
        threads = false,
        combine = :min,
    )[closure_firms]
    @test scores ≈ [
        0.87837245554794330,
        0.83928173885048851,
        0.84026525189212131,
        0.010491619358180508,
        0.83938875684430192,
        0.83919913726570627,
        0.83871134521356561,
        0.83862604469123081,
        0.83842687824673123,
        0.14782516685618610,
        0.013224577635286240,
        0.0064115243633385468,
        0.0046717737458951943,
        0.83930463691132073,
        0.011700530446422091,
        0.0086842594773570356,
    ] atol = 1e-12 rtol = 0
end

@testset "Internal helper functions" begin
    @test ESRIcascade._firm_selection(5, nothing) == [1, 2, 3, 4, 5]
    @test ESRIcascade._firm_selection(5, [5, 2, 3]) == [5, 2, 3]
    @test_throws BoundsError ESRIcascade._firm_selection(5, [0, 1])
    @test_throws ArgumentError ESRIcascade._firm_selection(5, [1, 1])
    @test ESRIcascade._linf_distance([1.0, 2.0, 3.0], [1.5, 1.0, 2.5]) == 1.0
end
