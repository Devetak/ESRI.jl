using Test
using Random
using SparseArrays
using LinearAlgebra: I

using ESRI

@testset "IndustryInfo" begin
    industry_ids = [1, 2, 1, 2]
    essential_industry = [true, false]
    info = IndustryInfo(industry_ids, essential_industry)

    @test ESRI.num_industries(info) == 2
    @test info.industry_of_firm == industry_ids
    @test info.essential_industry == essential_industry
    @test info.essential_firm == [true, false, true, false]

    @test_throws AssertionError IndustryInfo([0, 1], [true, false])
    @test_throws AssertionError IndustryInfo([3, 1], [true, false])
    @test_throws AssertionError IndustryInfo([1, 1], Bool[])
end

@testset "compute_esri (sparse) correctness + ranges" begin
    Random.seed!(123)
    N = 20
    W = sprand(N, N, 0.08) + 0.1I # ensure type stability + avoid empty graph

    industry_ids = rand(1:5, N)
    essential_industry = [true, true, true, false, false]
    info = IndustryInfo(industry_ids, essential_industry)

    esri_serial = compute_esri(W, info; maxiter = 30, tol = 1e-3, verbose = false, threads = false)
    esri_threaded = compute_esri(W, info; maxiter = 30, tol = 1e-3, verbose = false, threads = true)

    @test size(esri_serial) == (N,)
    @test all(isfinite, esri_serial)
    @test all(0 .<= esri_serial .<= 1 .+ 1e-6)

    # Results must be deterministic w.r.t. threading for this algorithm.
    @test esri_serial == esri_threaded
end

@testset "compute_esri (isolated firms edge case)" begin
    # Construct a sparse matrix where some firms have zero row-sum and/or column-sum.
    N = 8
    W = spzeros(N, N)
    W[1, 1] = 1.0
    W[2, 2] = 1.0
    W[3, 3] = 1.0
    W[4, 4] = 1.0
    W[1, 2] = 0.2
    W[2, 3] = 0.3

    industry_ids = [1, 2, 1, 2, 3, 3, 3, 3]
    essential_industry = [true, false, true]
    info = IndustryInfo(industry_ids, essential_industry)

    esri = compute_esri(W, info; maxiter = 20, tol = 1e-3, verbose = false, threads = false)
    @test all(isfinite, esri)
    @test all(0 .<= esri .<= 1 .+ 1e-6)
end

