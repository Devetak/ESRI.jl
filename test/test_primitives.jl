@testset "Impact matrix builders" begin
    W = [0.0 2.0; 3.0 0.0]
    upstream_dense = ESRIcascade.create_upstream_impact_matrix(W)
    upstream_sparse = ESRIcascade.create_upstream_impact_matrix(sparse(W))
    @test upstream_dense ≈ [0.0 1.0; 1.0 0.0] atol = 0 rtol = 0
    @test Matrix(upstream_sparse) ≈ upstream_dense atol = 0 rtol = 0

    W2 = [
        1.0 2.0 0.0
        0.0 4.0 0.0
        0.0 3.0 1.0
    ]
    info2 = IndustryInfo([1, 2, 1], [true, false])
    d1_dense, d2_dense = ESRIcascade.compute_downstream_impact_matrices(W2, info2)
    d1_sparse, d2_sparse = ESRIcascade.compute_downstream_impact_matrices(sparse(W2), info2)

    @test d1_dense[1, 2] ≈ 2 / 5
    @test d1_dense[3, 2] ≈ 3 / 5
    @test d1_dense[2, 2] == 0
    @test d2_dense[2, 2] ≈ 4 / 9
    @test d2_dense[1, 2] == 0
    @test d2_dense[3, 2] == 0
    @test Matrix(d1_sparse) ≈ d1_dense atol = 1e-12 rtol = 0
    @test Matrix(d2_sparse) ≈ d2_dense atol = 1e-12 rtol = 0

    W3 = sparse([1, 3], [2, 1], [2.0, 4.0], 4, 4)
    info3 = IndustryInfo([1, 2, 1, 2], [true, false])
    u3_dense = ESRIcascade.create_upstream_impact_matrix(Matrix(W3))
    u3_sparse = ESRIcascade.create_upstream_impact_matrix(W3)
    e3_dense, n3_dense = ESRIcascade.compute_downstream_impact_matrices(Matrix(W3), info3)
    e3_sparse, n3_sparse = ESRIcascade.compute_downstream_impact_matrices(W3, info3)

    @test Matrix(u3_sparse) ≈ u3_dense atol = 1e-12 rtol = 0
    @test Matrix(e3_sparse) ≈ e3_dense atol = 1e-12 rtol = 0
    @test Matrix(n3_sparse) ≈ n3_dense atol = 1e-12 rtol = 0
    @test Matrix(u3_sparse)[:, 2] == zeros(4)
    @test Matrix(e3_sparse)[:, 4] == zeros(4)
    @test Matrix(n3_sparse)[:, 4] == zeros(4)
end

@testset "Per-customer input classifications" begin
    W = [
        0.0 0.0 0.0 2.0
        0.0 0.0 0.0 3.0
        0.0 0.0 0.0 5.0
        0.0 0.0 0.0 7.0
    ]
    info = IndustryInfo(
        [1, 1, 2, 3],
        [
            0 0 2
            0 0 1
            1 0 0
        ],
    )

    essential_dense, nonessential_dense =
        ESRIcascade.compute_downstream_impact_matrices(W, info)
    essential_sparse, nonessential_sparse =
        ESRIcascade.compute_downstream_impact_matrices(sparse(W), info)

    @test essential_dense[1, 4] == 2 / 5
    @test essential_dense[2, 4] == 3 / 5
    @test essential_dense[3, 4] == 0
    @test essential_dense[4, 4] == 0
    @test nonessential_dense[3, 4] == 5 / 17
    @test nonessential_dense[1, 4] == 0
    @test nonessential_dense[2, 4] == 0
    @test nonessential_dense[4, 4] == 0
    @test Matrix(essential_sparse) == essential_dense
    @test Matrix(nonessential_sparse) == nonessential_dense
    @test length(essential_sparse.nzval) == 2
    @test length(nonessential_sparse.nzval) == 1
end

@testset "Legacy Boolean classification parity" begin
    W = [
        1.0 2.0 0.0 0.0
        0.0 1.0 3.0 0.0
        4.0 0.0 1.0 5.0
        0.0 6.0 0.0 1.0
    ]
    industries = [1, 1, 2, 3]
    legacy = IndustryInfo(industries, [true, false, true])
    classified = IndustryInfo(
        industries,
        [
            2 2 2
            1 1 1
            2 2 2
        ],
    )
    @test legacy.input_classification == classified.input_classification

    for weights in (W, sparse(W))
        legacy_essential, legacy_nonessential =
            ESRIcascade.compute_downstream_impact_matrices(weights, legacy)
        classified_essential, classified_nonessential =
            ESRIcascade.compute_downstream_impact_matrices(weights, classified)
        @test Matrix(classified_essential) == Matrix(legacy_essential)
        @test Matrix(classified_nonessential) == Matrix(legacy_nonessential)
    end

    legacy_econ = ESRIEconomy(sparse(W), legacy)
    classified_econ = ESRIEconomy(sparse(W), classified)
    @test esri(classified_econ; maxiter = 60, tol = 1e-10) ==
          esri(legacy_econ; maxiter = 60, tol = 1e-10)

    shock = [0.0, 1.0, 0.6, 1.0]
    legacy_result =
        esri_shock(legacy_econ, shock; details = true, maxiter = 60, tol = 1e-10)
    classified_result =
        esri_shock(classified_econ, shock; details = true, maxiter = 60, tol = 1e-10)
    @test classified_result.esri == legacy_result.esri
    @test classified_result.upstream == legacy_result.upstream
    @test classified_result.downstream == legacy_result.downstream

    default_econ = ESRIEconomy(sparse(W), IndustryInfo(industries))
    explicit_default_econ = ESRIEconomy(sparse(W), IndustryInfo(industries, fill(2, 3, 3)))
    @test esri(default_econ; maxiter = 60, tol = 1e-10) ==
          esri(explicit_default_econ; maxiter = 60, tol = 1e-10)
end

@testset "C++ mixed-classification parity" begin
    W = [
        0.0 3.0 2.0 7.0 5.0 11.0
        4.0 0.0 9.0 1.0 8.0 6.0
        10.0 2.0 0.0 4.0 3.0 7.0
        6.0 5.0 1.0 0.0 9.0 2.0
        8.0 3.0 6.0 10.0 0.0 4.0
        2.0 7.0 5.0 3.0 11.0 0.0
    ]
    info = IndustryInfo(
        [1, 2, 3, 1, 2, 3],
        [
            2 1 0
            0 2 1
            1 0 2
        ],
    )

    dense_econ = ESRIEconomy(W, info)
    econ = ESRIEconomy(sparse(W), info)
    permuted = ESRIcascade._permute_sparse_economy(econ, [2, 1, 3, 4, 5, 6])
    @test permuted.info.input_classification == info.input_classification

    @test esri(econ; combine = :downstream, maxiter = 150, tol = 1e-12) ≈
          esri(dense_econ; combine = :downstream, maxiter = 150, tol = 1e-12) atol = 1e-12 rtol =
        0

    result =
        esri(econ, 1; details = true, combine = :downstream, maxiter = 150, tol = 1e-12)
    @test result.esri ≈ 0.50687988214742541 atol = 1e-12 rtol = 0
    @test result.downstream ≈
          [0.0, 0.6, 0.7894068197483480, 0.0, 0.6111111111111110, 0.8786670560685980] atol =
        1e-12 rtol = 0
end

@testset "Propagation kernels" begin
    info = IndustryInfo([1, 1, 2], [true, false])
    sigmas = zeros(3)
    temp = zeros(2)
    ESRIcascade.compute_sigmas!(sigmas, [2.0, 1.0, 3.0], [1.0, 0.5, 0.25], info, temp)
    @test sigmas ≈ [0.8, 0.4, 1.0] atol = 1e-12 rtol = 0

    hd = [0.8, 0.3, 1.0]
    ess = [1.0 0.0 0.0; 0.5 0.0 0.0; 0.0 0.0 0.0]
    non = [0.0 0.2 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0]
    emat = zeros(3, 2)
    nvec = zeros(3)
    ESRIcascade._accumulate_downstream_components!(
        emat,
        nvec,
        hd,
        [0.5, 1.0, 0.2],
        ess,
        non,
        info,
    )
    @test emat ≈ [0.45 0.0; 0.0 0.0; 0.0 0.0] atol = 1e-12 rtol = 0
    @test nvec ≈ [0.0, 0.02, 0.7] atol = 1e-12 rtol = 0

    out = zeros(3)
    ESRIcascade.downstream_step!(out, emat, nvec, [1.0, 0.5, 0.9])
    @test out ≈ [0.55, 0.5, 0.3] atol = 1e-12 rtol = 0

    curr_u = zeros(2)
    ESRIcascade.upstream_step!(
        curr_u,
        [0.0 0.5; 1.0 0.0],
        [0.6, 0.8],
        [1.0, 0.7],
        [1.0, 0.0],
    )
    @test curr_u ≈ [0.8, 0.7] atol = 1e-12 rtol = 0

    curr_u_sp = zeros(2)
    ESRIcascade.upstream_step!(
        curr_u_sp,
        sparse([0.0 0.5; 1.0 0.0]),
        [0.6, 0.8],
        [1.0, 0.7],
        [1.0, 0.0],
    )
    @test curr_u_sp ≈ curr_u atol = 1e-12 rtol = 0
end

@testset "Sparse downstream parity" begin
    W = [
        0.0 2.0 3.0 1.0
        4.0 0.0 5.0 2.0
        6.0 1.0 0.0 3.0
        2.0 4.0 1.0 0.0
    ]
    info = IndustryInfo([1, 1, 2, 2], [2 1; 1 2])
    psi = [0.0, 1.0, 0.4, 0.8]
    dense = ESRIEconomy(W, info)
    sparse_econ = ESRIEconomy(sparse(W), info)

    dense_scores = esri(dense; firm_indices = 1:4, maxiter = 80, tol = 1e-12)
    sparse_scores = esri(sparse_econ; firm_indices = 1:4, maxiter = 80, tol = 1e-12)
    @test sparse_scores ≈ dense_scores atol = 1e-12 rtol = 0

    dense_result = esri_shock(dense, psi; details = true, maxiter = 80, tol = 1e-12)
    sparse_result = esri_shock(sparse_econ, psi; details = true, maxiter = 80, tol = 1e-12)
    @test sparse_result.esri ≈ dense_result.esri atol = 1e-12 rtol = 0
    @test sparse_result.downstream ≈ dense_result.downstream atol = 1e-12 rtol = 0
end

@testset "Integer helper inputs" begin
    W = [1 1; 0 2]
    info = IndustryInfo([1, 1], [true])

    upstream = ESRIcascade.create_upstream_impact_matrix(W)
    essential, nonessential = ESRIcascade.compute_downstream_impact_matrices(W, info)

    @test upstream ≈ [0.5 0.0; 0.5 1.0] atol = 1e-12 rtol = 0
    @test essential ≈ [1.0 1 / 3; 0.0 2 / 3] atol = 1e-12 rtol = 0
    @test nonessential ≈ zeros(2, 2) atol = 1e-12 rtol = 0
end
