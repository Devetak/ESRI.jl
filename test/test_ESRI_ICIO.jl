using ESRIcascade
using DelimitedFiles
using SparseArrays
using Test

data_dir = normpath(joinpath(@__DIR__, "..", "data", "ICIO"))
codes = ihs_industry_codes()
production_codes = vec(readdlm(joinpath(data_dir, "p_icio.csv"), ',', String))
weight_matrix = sparse(readdlm(joinpath(data_dir, "W_ICIO.csv"), ',', Float64))
reference = readdlm(joinpath(data_dir, "ESRI_icio.csv"), ',', Float64)

nace2(code) = code[1:2]
industry_by_nace2 = Dict{String,Int}()
for (industry, code) in pairs(codes)
    get!(industry_by_nace2, nace2(code), industry)
end
industry_ids = [get(industry_by_nace2, nace2(code), 0) for code in production_codes]

@testset "ICIO reference data" begin
    @test size(weight_matrix) == (length(production_codes), length(production_codes))
    @test size(reference) == (length(production_codes), 3)
    @test all(industry_ids .> 0)
    @test all(isfinite, reference)
    @test all(reference .>= 0)

    economy = ESRIEconomy(weight_matrix, IndustryInfo(industry_ids, ihs_input_classification()))
    for (combine, column) in ((:min, 1), (:downstream, 2), (:upstream, 3))
        scores = esri(economy; combine = combine)
        @test isapprox(scores, reference[:, column]; atol = 1e-6, rtol = 0)
    end
end
