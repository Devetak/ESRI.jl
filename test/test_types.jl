@testset "IndustryInfo" begin
    industry_ids = [1, 2, 1, 2]
    essential_industry = [true, false]
    info = IndustryInfo(industry_ids, essential_industry)

    @test ESRI.num_industries(info) == 2
    @test info.industry_of_firm == industry_ids
    @test info.essential_industry == essential_industry
    @test info.essential_firm == [true, false, true, false]

    @test IndustryInfo(Int32[1, 2, 1, 2], Bool[true, false]).industry_of_firm == [1, 2, 1, 2]
    @test_throws ArgumentError IndustryInfo([0, 1], [true, false])
    @test_throws ArgumentError IndustryInfo([3, 1], [true, false])
    @test_throws ArgumentError IndustryInfo([1, 1], Bool[])
end

@testset "Internal helper functions" begin
    @test ESRI._firm_selection(5, nothing) == [1, 2, 3, 4, 5]
    @test ESRI._firm_selection(5, [5, 2, 3]) == [5, 2, 3]
    @test_throws BoundsError ESRI._firm_selection(5, [0, 1])
    @test_throws ArgumentError ESRI._firm_selection(5, [1, 1])
    @test ESRI._linf_distance([1.0, 2.0, 3.0], [1.5, 1.0, 2.5]) == 1.0
end
