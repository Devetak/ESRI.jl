@testset "Economy-wide API invariants" begin
    W_sparse, info = sparse_fixture()
    W_dense = Matrix(W_sparse)
    econ = ESRIEconomy(W_sparse, info)

    old_blas = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        scores_sparse = esri(econ; maxiter = 30, tol = 1e-3, threads = false)
        scores_dense = compute_esri(W_dense, info; maxiter = 30, tol = 1e-3, threads = false)
        scores_threaded = esri(econ; maxiter = 30, tol = 1e-3, threads = true)

        @test length(econ) == size(W_sparse, 1) == econ.n
        @test size(econ.upstream_impact) == size(W_sparse)
        @test size(econ.downstream_impact_essential) == size(W_sparse)
        @test size(econ.downstream_impact_nonessential) == size(W_sparse)
        assert_bounded(scores_sparse; atol = 1e-6)
        assert_bounded(scores_dense; atol = 1e-6)
        @test scores_sparse ≈ scores_dense atol = 1e-6 rtol = 1e-6
        @test scores_sparse == scores_threaded

        subset = [2, 5, 9]
        subset_scores = esri(econ; maxiter = 30, tol = 1e-3, firm_indices = subset)
        @test all(subset_scores[setdiff(1:econ.n, subset)] .== 0)
        @test subset_scores[subset] ≈ scores_sparse[subset] atol = 1e-12 rtol = 1e-12
        @test all(iszero, esri(econ; maxiter = 30, tol = 1e-3, firm_indices = Int[]))

        @test_throws ArgumentError esri(econ; firm_indices = [1, 1])
    finally
        BLAS.set_num_threads(old_blas)
    end
end

@testset "Single-scenario API invariants" begin
    W, info = deterministic_fixture()
    econ = ESRIEconomy(W, info)

    value = esri(econ, 1; maxiter = 100, tol = 1e-12)
    details = esri(econ, 1; details = true, maxiter = 100, tol = 1e-12)
    up_only = esri(econ, 1; components = :upstream, maxiter = 100, tol = 1e-12)
    down_only = esri(econ, 1; components = :downstream, maxiter = 100, tol = 1e-12)

    @test details isa ESRIResult
    @test details.esri ≈ value atol = 1e-12 rtol = 0
    @test up_only.esri ≈ value atol = 1e-12 rtol = 0
    @test down_only.esri ≈ value atol = 1e-12 rtol = 0
    @test up_only.upstream == details.upstream
    @test down_only.downstream == details.downstream
    assert_bounded(details.upstream)
    assert_bounded(details.downstream)

    psi = [0.2, 0.8, 1.0]
    direct = esri(econ, 1; shock = psi, details = true, maxiter = 100, tol = 1e-12)
    wrapped = esri_shock(econ, psi; details = true, maxiter = 100, tol = 1e-12)
    @test direct.esri ≈ wrapped.esri atol = 1e-12 rtol = 0
    @test direct.upstream ≈ wrapped.upstream atol = 1e-12 rtol = 0
    @test direct.downstream ≈ wrapped.downstream atol = 1e-12 rtol = 0
    @test all(direct.upstream .<= psi .+ 1e-12)
    @test all(direct.downstream .<= psi .+ 1e-12)

    W_sparse, info_sparse = sparse_fixture(seed = 91, n = 14, density = 0.1, nindustries = 3)
    econ_sparse = ESRIEconomy(W_sparse, info_sparse)
    psi_sparse = fill(1.0, 14)
    psi_sparse[[1, 4, 7]] .= (0.1, 0.6, 0.0)
    direct_sparse = esri(econ_sparse, 4; shock = psi_sparse, details = true, maxiter = 25, tol = 1e-3)
    wrapped_sparse = esri_shock(econ_sparse, psi_sparse; details = true, maxiter = 25, tol = 1e-3)
    @test direct_sparse.esri ≈ wrapped_sparse.esri atol = 1e-12 rtol = 0
    @test direct_sparse.upstream ≈ wrapped_sparse.upstream atol = 1e-12 rtol = 0
    @test direct_sparse.downstream ≈ wrapped_sparse.downstream atol = 1e-12 rtol = 0
end

@testset "Hand-computed and zero-output cases" begin
    W = [
        1.0 1.0
        0.0 1.0
    ]
    econ = ESRIEconomy(W, IndustryInfo([1, 1], [true]))
    details = esri(econ, 1; details = true, maxiter = 100, tol = 1e-12)

    @test details.upstream ≈ [0.0, 1.0] atol = 1e-12 rtol = 0
    @test details.downstream ≈ [0.0, 0.0] atol = 1e-10 rtol = 0
    @test esri(econ, 1; combine = :upstream, maxiter = 100, tol = 1e-12) ≈ 2 / 3 atol = 1e-12 rtol = 0
    @test esri(econ, 1; combine = :downstream, maxiter = 100, tol = 1e-12) ≈ 1.0 atol = 1e-10 rtol = 0
    @test details.esri ≈ 1.0 atol = 1e-10 rtol = 0

    econ0 = ESRIEconomy(zeros(4, 4), IndustryInfo([1, 2, 1, 2], [true, false]))
    @test econ0.total_output == 0.0
    @test esri(econ0; maxiter = 20, tol = 1e-6) == zeros(4)
    d = esri(econ0, 2; details = true, maxiter = 20, tol = 1e-6)
    @test d.esri == 0.0
    @test d.upstream == [1.0, 0.0, 1.0, 1.0]
    @test d.downstream == [1.0, 0.0, 1.0, 1.0]
    @test esri(econ0, 2; final_weights = ones(4), maxiter = 20, tol = 1e-6) == 1.0

    econ0_sparse = ESRIEconomy(spzeros(4, 4), IndustryInfo([1, 2, 1, 2], [true, false]))
    @test econ0_sparse.total_output == 0.0
    @test esri(econ0_sparse; maxiter = 20, tol = 1e-6) == zeros(4)
    d_sparse = esri(econ0_sparse, 2; details = true, maxiter = 20, tol = 1e-6)
    @test d_sparse.esri == 0.0
    @test d_sparse.upstream == [1.0, 0.0, 1.0, 1.0]
    @test d_sparse.downstream == [1.0, 0.0, 1.0, 1.0]
    @test esri(econ0_sparse, 2; final_weights = ones(4), maxiter = 20, tol = 1e-6) == 1.0
end

@testset "Validation and logging" begin
    W = [1.0 0.1; 0.0 1.0]
    info = IndustryInfo([1, 1], [true])
    econ = ESRIEconomy(W, info)

    @test_throws BoundsError esri(econ, 0)
    @test_throws ArgumentError esri(econ; combine = :bad)
    @test_throws ArgumentError esri(econ, 1; components = :bad)
    @test_throws DimensionMismatch esri(econ; final_weights = [1.0])
    @test_throws DimensionMismatch esri(econ, 1; shock = [0.0])
    @test_throws DimensionMismatch ESRIEconomy(ones(2, 3), info)
    @test_throws DomainError esri(econ, 1; shock = [NaN, 1.0])
    @test_throws DomainError esri_shock(econ, [1.0, Inf])
    @test_throws DomainError esri(econ; final_weights = [1.0, -0.2])
    @test_throws DomainError ESRIEconomy([1.0 NaN; 0.0 1.0], info)

    @test_logs (:info, r"joint iteration") begin
        @test isfinite(esri(econ, 1; verbose = true, maxiter = 10, tol = 0.0))
    end

    ref = compute_esri(W, info; maxiter = 10, tol = 1e-2, threads = true)
    if Threads.nthreads() > 1
        @test_logs (:warn, r"Ignoring `verbose=true` because progress UI is disabled in threaded mode\.") begin
            @test compute_esri(W, info; maxiter = 10, tol = 1e-2, verbose = true, threads = true) == ref
        end
    else
        @test compute_esri(W, info; maxiter = 10, tol = 1e-2, verbose = true, threads = true) == ref
    end
end
