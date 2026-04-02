@testset "compute_esri (sparse) correctness + ranges" begin
    Random.seed!(123)
    n = 20
    W = sprand(n, n, 0.08) + 0.1I

    industry_ids = rand(1:5, n)
    essential_industry = [true, true, true, false, false]
    info = IndustryInfo(industry_ids, essential_industry)

    esri_serial = compute_esri(W, info; maxiter = 30, tol = 1e-3, verbose = false, threads = false)
    esri_threaded = compute_esri(W, info; maxiter = 30, tol = 1e-3, verbose = false, threads = true)

    @test size(esri_serial) == (n,)
    @test all(isfinite, esri_serial)
    @test all(0 .<= esri_serial .<= 1 .+ 1e-6)
    @test esri_serial == esri_threaded
end

@testset "compute_esri (isolated firms edge case)" begin
    n = 8
    W = spzeros(n, n)
    W[1, 1] = 1.0
    W[2, 2] = 1.0
    W[3, 3] = 1.0
    W[4, 4] = 1.0
    W[1, 2] = 0.2
    W[2, 3] = 0.3

    industry_ids = [1, 2, 1, 2, 3, 3, 3, 3]
    essential_industry = [true, false, true]
    info = IndustryInfo(industry_ids, essential_industry)

    values = compute_esri(W, info; maxiter = 20, tol = 1e-3, verbose = false, threads = false)
    @test all(isfinite, values)
    @test all(0 .<= values .<= 1 .+ 1e-6)
end

@testset "compute_esri sparse vs dense match" begin
    old_n = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        Random.seed!(202)
        n = 25
        W_sparse = sprand(n, n, 0.12) + 0.1I
        W_dense = Matrix(W_sparse)

        industry_ids = rand(1:4, n)
        essential_industry = [true, false, true, false]
        info = IndustryInfo(industry_ids, essential_industry)

        esri_sparse = compute_esri(W_sparse, info; maxiter = 25, tol = 1e-3, verbose = false, threads = false)
        esri_dense = compute_esri(W_dense, info; maxiter = 25, tol = 1e-3, verbose = false, threads = false)

        @test size(esri_sparse) == (n,)
        @test size(esri_dense) == (n,)
        @test all(0 .<= esri_sparse .<= 1 .+ 1e-6)
        @test all(0 .<= esri_dense .<= 1 .+ 1e-6)
        @test esri_dense ≈ esri_sparse atol = 1e-6 rtol = 1e-6
    finally
        BLAS.set_num_threads(old_n)
    end
end

@testset "compute_esri (firm_indices subset + validation)" begin
    n = 12
    W = spzeros(n, n)
    @inbounds for i in 1:n
        W[i, i] = 1.0
    end
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
    essential_industry = [true, false, true]
    info = IndustryInfo(industry_ids, essential_industry)

    subset = [2, 5, 9]
    full = compute_esri(W, info; maxiter = 25, tol = 1e-3, verbose = false, threads = false)
    partial = compute_esri(W, info; maxiter = 25, tol = 1e-3, verbose = false, threads = false, firm_indices = subset)

    @test all(partial[setdiff(1:n, subset)] .== 0)
    @test partial[subset] == full[subset]

    partial_threaded = compute_esri(W, info; maxiter = 25, tol = 1e-3, verbose = false, threads = true, firm_indices = subset)
    @test partial_threaded[subset] == full[subset]

    @test_throws ArgumentError compute_esri(W, info; maxiter = 5, tol = 1e-2, verbose = false, threads = false, firm_indices = [1, 1])
end

@testset "compute_esri (verbose in threaded mode logs)" begin
    Random.seed!(808)
    n = 6
    W = sprand(n, n, 0.3) + 0.1I
    industry_ids = rand(1:2, n)
    essential_industry = [true, false]
    info = IndustryInfo(industry_ids, essential_industry)

    ref = compute_esri(W, info; maxiter = 10, tol = 1e-2, verbose = false, threads = true)
    if Threads.nthreads() > 1
        @test_logs (:warn, r"Ignoring `verbose=true` because progress UI is disabled in threaded mode\.") begin
            got = compute_esri(W, info; maxiter = 10, tol = 1e-2, verbose = true, threads = true)
            @test got == ref
        end
    else
        got = compute_esri(W, info; maxiter = 10, tol = 1e-2, verbose = true, threads = true)
        @test got == ref
    end
end

@testset "ESRIEconomy + esri API" begin
    Random.seed!(404)
    n = 18
    W = sprand(n, n, 0.12) + 0.1I
    industry_ids = rand(1:4, n)
    essential_industry = [true, false, true, false]
    info = IndustryInfo(industry_ids, essential_industry)

    econ = ESRIEconomy(W, info)
    @test length(econ) == n
    @test econ.n == n
    @test size(econ.upstream_impact) == (n, n)
    @test size(econ.downstream_impact_essential) == (n, n)
    @test size(econ.downstream_impact_nonessential) == (n, n)
    @test length(econ.column_sums) == n
    @test length(econ.row_sums) == n

    old_api = compute_esri(W, info; maxiter = 25, tol = 1e-3, verbose = false, threads = false)
    new_api = esri(econ; maxiter = 25, tol = 1e-3, verbose = false, threads = false)
    @test new_api ≈ old_api atol = 1e-10 rtol = 1e-10

    subset = [2, 7, 11]
    subset_values = esri(econ; maxiter = 25, tol = 1e-3, verbose = false, threads = false, firm_indices = subset)
    @test all(subset_values[setdiff(1:n, subset)] .== 0)
    @test subset_values[subset] == new_api[subset]

    i = 7
    single = esri(econ, i; maxiter = 25, tol = 1e-3, verbose = false)
    @test single ≈ new_api[i] atol = 1e-10 rtol = 1e-10

    details = esri(econ, i; maxiter = 25, tol = 1e-3, verbose = false, details = true)
    @test details isa ESRIResult
    @test details.esri ≈ new_api[i] atol = 1e-10 rtol = 1e-10
    @test length(details.upstream) == n
    @test length(details.downstream) == n
    @test all(0 .<= details.upstream .<= 1 .+ 1e-8)
    @test all(0 .<= details.downstream .<= 1 .+ 1e-8)

    only_upstream = esri(econ, i; maxiter = 25, tol = 1e-3, components = :upstream)
    only_downstream = esri(econ, i; maxiter = 25, tol = 1e-3, components = :downstream)
    @test haskey(only_upstream, :upstream)
    @test !haskey(only_upstream, :downstream)
    @test haskey(only_downstream, :downstream)
    @test !haskey(only_downstream, :upstream)
end

@testset "public API validation errors" begin
    W = [1.0 0.1; 0.0 1.0]
    info = IndustryInfo([1, 1], [true])
    econ = ESRIEconomy(W, info)

    @test_throws BoundsError esri(econ, 0)
    @test_throws ArgumentError esri(econ; combine = :bad)
    @test_throws ArgumentError esri(econ, 1; components = :bad)
    @test_throws DimensionMismatch esri(econ; final_weights = [1.0])
    @test_throws DimensionMismatch esri(econ, 1; shock = [0.0])
    @test_throws DomainError esri(econ, 1; shock = [1.2, 0.3])
    @test_throws DimensionMismatch ESRIEconomy(ones(2, 3), info)
end
