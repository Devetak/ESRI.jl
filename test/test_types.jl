@testset "IndustryInfo" begin
    industry_ids = [1, 2, 1, 2]
    essential_industry = [true, false]
    info = IndustryInfo(industry_ids, essential_industry)

    @test ESRIcascade.num_industries(info) == 2
    @test info.industry_of_firm == industry_ids
    @test info.essential_industry == essential_industry
    @test info.essential_firm == [true, false, true, false]

    default_info = IndustryInfo(industry_ids)
    @test default_info.input_classification == fill(UInt8(2), 2, 2)
    @test default_info.essential_firm == trues(length(industry_ids))
    @test_throws ArgumentError IndustryInfo(Int[])

    @test IndustryInfo(Int32[1, 2, 1, 2], Bool[true, false]).industry_of_firm == [1, 2, 1, 2]
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
end

@testset "Bundled IHS classification" begin
    classification = ihs_input_classification()
    @test size(classification) == (616, 616)
    @test ihs_industry_codes()[[1, end]] == ["0111", "9999"]
    @test count(==(UInt8(0)), classification) == 269204
    @test count(==(UInt8(1)), classification) == 45013
    @test count(==(UInt8(2)), classification) == 65239

    classification[1, 1] = 0
    @test ihs_input_classification()[1, 1] == 2
end

@testset "Internal helper functions" begin
    @test ESRIcascade._firm_selection(5, nothing) == [1, 2, 3, 4, 5]
    @test ESRIcascade._firm_selection(5, [5, 2, 3]) == [5, 2, 3]
    @test_throws BoundsError ESRIcascade._firm_selection(5, [0, 1])
    @test_throws ArgumentError ESRIcascade._firm_selection(5, [1, 1])
    @test ESRIcascade._linf_distance([1.0, 2.0, 3.0], [1.5, 1.0, 2.5]) == 1.0
end
