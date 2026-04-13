@testset "Custom weights and combine modes" begin
    W, info = deterministic_fixture()
    econ = ESRIEconomy(W, info)
    psi = [0.2, 0.8, 1.0]
    weights = [5.0, 7.0, 11.0]

    scenario = esri_shock(econ, psi; details = true, maxiter = 100, tol = 1e-12)
    up_only = esri_shock(econ, psi; components = :upstream, combine = :upstream, maxiter = 100, tol = 1e-12)
    down_only = esri_shock(econ, psi; components = :downstream, combine = :downstream, maxiter = 100, tol = 1e-12)

    expected_up = sum(econ.row_sums .* (1 .- up_only.upstream)) / econ.total_output
    expected_down = sum(econ.row_sums .* (1 .- down_only.downstream)) / econ.total_output

    @test up_only.esri ≈ expected_up atol = 1e-12 rtol = 0
    @test down_only.esri ≈ expected_down atol = 1e-12 rtol = 0

    scalar = esri_shock(econ, psi; final_weights = weights, maxiter = 100, tol = 1e-12)
    detail = esri_shock(econ, psi; final_weights = weights, details = true, maxiter = 100, tol = 1e-12)
    expected = sum(weights .* (1 .- min.(detail.upstream, detail.downstream))) / econ.total_output

    @test scalar ≈ expected atol = 1e-12 rtol = 0
    @test detail.esri ≈ expected atol = 1e-12 rtol = 0
    @test compute_esri_shock(W, info, psi; final_weights = weights, maxiter = 100, tol = 1e-12) ≈ expected atol = 1e-12 rtol = 0
    @test compute_esri(W, info, 1; shock = psi, final_weights = weights, maxiter = 100, tol = 1e-12) ≈ expected atol = 1e-12 rtol = 0
end

@testset "Dense, sparse, and threaded option parity" begin
    W_sparse = sparse(
        [1, 1, 2, 2, 3, 4, 5, 6, 3, 4],
        [1, 2, 2, 3, 3, 4, 5, 6, 4, 5],
        [1.0, 0.4, 1.0, 0.2, 1.0, 1.0, 1.0, 1.0, 0.3, 0.1],
        6,
        6,
    )
    W_dense = Matrix(W_sparse)
    info = IndustryInfo([1, 1, 2, 2, 1, 2], [true, false])
    weights = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]

    sparse_scores = compute_esri(W_sparse, info; final_weights = weights, combine = :downstream, maxiter = 80, tol = 1e-10, threads = false)
    dense_scores = compute_esri(W_dense, info; final_weights = weights, combine = :downstream, maxiter = 80, tol = 1e-10, threads = false)
    threaded_scores = compute_esri(W_sparse, info; final_weights = weights, combine = :downstream, maxiter = 80, tol = 1e-10, threads = true)

    @test sparse_scores ≈ dense_scores atol = 1e-12 rtol = 0
    @test sparse_scores ≈ threaded_scores atol = 1e-12 rtol = 0
end
