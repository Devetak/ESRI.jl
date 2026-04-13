@testset "Coupled psi cap and upstream boundary" begin
    W = [
        1.0 0.0 0.0
        1.0 1.0 0.0
        0.0 0.0 0.0
    ]
    info = IndustryInfo([1, 1, 2], [true, false])
    econ = ESRIEconomy(W, info)

    details = esri(econ, 1; maxiter = 50, tol = 1e-6, details = true)
    @test details.downstream[1] == 0.0
    @test details.upstream[1] == 0.0
    @test details.upstream[3] == 1.0
    @test all(0 .<= details.downstream .<= 1 .+ 1e-8)
    @test all(0 .<= details.upstream .<= 1 .+ 1e-8)
end

@testset "Wrapper consistency for custom shock" begin
    Random.seed!(91)
    n = 14
    W = sprand(n, n, 0.1) + 0.1I
    info = IndustryInfo(rand(1:3, n), [true, false, true])
    econ = ESRIEconomy(W, info)

    psi = rand(n)
    direct = esri_shock(econ, psi; maxiter = 25, tol = 1e-3)
    wrapped = compute_esri_shock(W, info, psi; maxiter = 25, tol = 1e-3)
    @test direct ≈ wrapped atol = 1e-12 rtol = 0
end

@testset "Split convergence matches lockstep reference" begin
    W_dense, info_dense = deterministic_fixture()
    econ_dense = ESRIEconomy(W_dense, info_dense)

    ref_dense = reference_lockstep_scenario(econ_dense, 1; maxiter = 100, tol = 1e-12)
    got_dense = esri(econ_dense, 1; details = true, maxiter = 100, tol = 1e-12)
    @test got_dense.esri ≈ ref_dense.esri atol = 1e-12 rtol = 0
    @test got_dense.upstream ≈ ref_dense.upstream atol = 1e-12 rtol = 0
    @test got_dense.downstream ≈ ref_dense.downstream atol = 1e-12 rtol = 0

    W_sparse, info_sparse, _ = powerlaw_fixture(seed = 8, n = 40, mean_degree = 6, max_degree = 32)
    econ_sparse = ESRIEconomy(W_sparse, info_sparse)
    psi = fill(1.0, 40)
    psi[[2, 9, 17]] .= (0.0, 0.4, 0.8)

    ref_sparse = reference_lockstep_scenario(econ_sparse, 9; maxiter = 60, tol = 1e-6, shock = psi)
    got_sparse = esri(econ_sparse, 9; details = true, maxiter = 60, tol = 1e-6, shock = psi)
    @test got_sparse.esri ≈ ref_sparse.esri atol = 1e-10 rtol = 0
    @test got_sparse.upstream ≈ ref_sparse.upstream atol = 1e-10 rtol = 0
    @test got_sparse.downstream ≈ ref_sparse.downstream atol = 1e-6 rtol = 1e-6
end

@testset "Power-law sparse full-economy parity" begin
    W_sparse, info, _ = powerlaw_fixture(seed = 17, n = 60, mean_degree = 7, max_degree = 48, nindustries = 8)
    W_dense = Matrix(W_sparse)
    econ = ESRIEconomy(W_sparse, info)

    sparse_scores = esri(econ; maxiter = 35, tol = 1e-3, threads = false)
    explicit_scores = esri(econ; maxiter = 35, tol = 1e-3, threads = false, firm_indices = collect(1:size(W_sparse, 1)))
    dense_scores = compute_esri(W_dense, info; maxiter = 35, tol = 1e-3, threads = false)
    threaded_scores = compute_esri(W_sparse, info; maxiter = 35, tol = 1e-3, threads = true)

    @test sparse_scores ≈ explicit_scores atol = 1e-10 rtol = 1e-10
    @test sparse_scores ≈ dense_scores atol = 1e-8 rtol = 1e-8
    @test sparse_scores == threaded_scores
end
