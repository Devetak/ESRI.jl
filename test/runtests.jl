using Test
using Random
using SparseArrays
using LinearAlgebra: I, BLAS

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

@testset "compute_esri sparse vs dense match" begin
    # Keep BLAS single-threaded to reduce run-to-run numeric noise.
    old_n = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        Random.seed!(202)
        N = 25
        W_sparse = sprand(N, N, 0.12) + 0.1I
        W_dense = Matrix(W_sparse)

        industry_ids = rand(1:4, N)
        essential_industry = [true, false, true, false]
        info = IndustryInfo(industry_ids, essential_industry)

        esri_sparse = compute_esri(W_sparse, info; maxiter = 25, tol = 1e-3, verbose = false, threads = false)
        esri_dense = compute_esri(W_dense, info; maxiter = 25, tol = 1e-3, verbose = false, threads = false)

        @test esri_sparse isa AbstractVector
        @test esri_dense isa AbstractVector
        @test size(esri_sparse) == (N,)
        @test size(esri_dense) == (N,)
        @test all(0 .<= esri_sparse .<= 1 .+ 1e-6)
        @test all(0 .<= esri_dense .<= 1 .+ 1e-6)

        # Dense path uses BLAS `mul!` while sparse path uses explicit sparse loops.
        # Results should be numerically very close for this deterministic fixture.
        @test esri_dense ≈ esri_sparse atol = 1e-6 rtol = 1e-6
    finally
        BLAS.set_num_threads(old_n)
    end
end

@testset "compute_esri (firm_indices subset + validation)" begin
    N = 12
    W = spzeros(N, N)
    # Deterministic sparse directed network.
    W[1, 1] = 1.0
    W[2, 2] = 1.0
    W[3, 3] = 1.0
    W[4, 4] = 1.0
    W[5, 5] = 1.0
    W[6, 6] = 1.0
    W[7, 7] = 1.0
    W[8, 8] = 1.0
    W[9, 9] = 1.0
    W[10, 10] = 1.0
    W[11, 11] = 1.0
    W[12, 12] = 1.0

    # Cross edges (kept sparse)
    W[1, 2] = 0.2
    W[2, 3] = 0.3
    W[3, 5] = 0.4
    W[4, 6] = 0.1
    W[5, 7] = 0.25
    W[6, 8] = 0.15
    W[7, 9] = 0.05
    W[8, 10] = 0.35
    W[9, 11] = 0.2
    W[10, 12] = 0.1

    industry_ids = [1, 2, 1, 2, 3, 3, 2, 1, 3, 2, 1, 2]
    essential_industry = [true, false, true] # 3 industries
    info = IndustryInfo(industry_ids, essential_industry)

    subset = [2, 5, 9]
    esri_full = compute_esri(W, info; maxiter = 25, tol = 1e-3, verbose = false, threads = false)
    esri_sub = compute_esri(W, info; maxiter = 25, tol = 1e-3, verbose = false, threads = false, firm_indices = subset)

    # Uncomputed entries must remain zero.
    remaining = setdiff(1:N, subset)
    @test all(esri_sub[remaining] .== 0)
    # Computed entries must match the full run.
    @test esri_sub[subset] == esri_full[subset]

    # Threaded mode must preserve the same results.
    esri_sub_threaded = compute_esri(W, info; maxiter = 25, tol = 1e-3, verbose = false, threads = true, firm_indices = subset)
    @test esri_sub_threaded[subset] == esri_full[subset]

    # `firm_indices` must be unique.
    @test_throws AssertionError compute_esri(W, info; maxiter = 5, tol = 1e-2, verbose = false, threads = false, firm_indices = [1, 1])
end

@testset "compute_esri (verbose in threaded mode logs)" begin
    # The threaded path disables progress UI; it should warn (not error).
    N = 6
    W = sprand(N, N, 0.3) + 0.1I
    industry_ids = rand(1:2, N)
    essential_industry = [true, false]
    info = IndustryInfo(industry_ids, essential_industry)

    esri_serial = compute_esri(W, info; maxiter = 10, tol = 1e-2, verbose = false, threads = false)
    esri_threaded_ref = compute_esri(W, info; maxiter = 10, tol = 1e-2, verbose = false, threads = true)

    @test_logs (:warn, r"Ignoring `verbose=true` because progress UI is disabled in threaded mode\." ) begin
        esri_threaded = compute_esri(W, info; maxiter = 10, tol = 1e-2, verbose = true, threads = true)
        @test esri_threaded == esri_threaded_ref
    end

    @test esri_threaded_ref == compute_esri(W, info; maxiter = 10, tol = 1e-2, verbose = false, threads = false)
    @test all(0 .<= esri_serial .<= 1 .+ 1e-6)
end

