@testset "Impact matrix builders" begin
    W = [0.0 2.0; 3.0 0.0]
    upstream_dense = ESRI.create_upstream_impact_matrix(W)
    upstream_sparse = ESRI.create_upstream_impact_matrix(sparse(W))
    @test upstream_dense ≈ [0.0 1.0; 1.0 0.0] atol = 0 rtol = 0
    @test Matrix(upstream_sparse) ≈ upstream_dense atol = 0 rtol = 0

    W2 = [
        1.0 2.0 0.0
        0.0 4.0 0.0
        0.0 3.0 1.0
    ]
    info2 = IndustryInfo([1, 2, 1], [true, false])
    d1_dense, d2_dense = ESRI.compute_downstream_impact_matrices(W2, info2)
    d1_sparse, d2_sparse = ESRI.compute_downstream_impact_matrices(sparse(W2), info2)

    @test d1_dense[1, 2] ≈ 2 / 5
    @test d1_dense[3, 2] ≈ 3 / 5
    @test d1_dense[2, 2] == 0
    @test d2_dense[2, 2] ≈ 4 / 9
    @test d2_dense[1, 2] == 0
    @test d2_dense[3, 2] == 0
    @test Matrix(d1_sparse) ≈ d1_dense atol = 1e-12 rtol = 0
    @test Matrix(d2_sparse) ≈ d2_dense atol = 1e-12 rtol = 0
end

@testset "Propagation kernels" begin
    info = IndustryInfo([1, 1, 2], [true, false])
    sigmas = zeros(3)
    temp = zeros(2)
    ESRI.compute_sigmas!(sigmas, [2.0, 1.0, 3.0], [1.0, 0.5, 0.25], info, temp)
    @test sigmas ≈ [0.8, 0.4, 1.0] atol = 1e-12 rtol = 0

    hd = [0.8, 0.3, 1.0]
    ess = [1.0 0.0 0.0; 0.5 0.0 0.0; 0.0 0.0 0.0]
    non = [0.0 0.2 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0]
    emat = zeros(3, 2)
    nvec = zeros(3)
    ESRI._accumulate_downstream_components!(emat, nvec, hd, [0.5, 1.0, 0.2], ess, non, info)
    @test emat ≈ [0.45 0.0; 0.0 0.0; 0.0 0.0] atol = 1e-12 rtol = 0
    @test nvec ≈ [0.0, 0.02, 0.7] atol = 1e-12 rtol = 0

    out = zeros(3)
    ESRI.downstream_step!(out, emat, nvec, [1.0, 0.5, 0.9])
    @test out ≈ [0.55, 0.5, 0.3] atol = 1e-12 rtol = 0

    curr_u = zeros(2)
    ESRI.upstream_step!(curr_u, [0.0 0.5; 1.0 0.0], [0.6, 0.8], [1.0, 0.7], [1.0, 0.0])
    @test curr_u ≈ [0.8, 0.7] atol = 1e-12 rtol = 0

    curr_u_sp = zeros(2)
    ESRI.upstream_step!(curr_u_sp, sparse([0.0 0.5; 1.0 0.0]), [0.6, 0.8], [1.0, 0.7], [1.0, 0.0])
    @test curr_u_sp ≈ curr_u atol = 1e-12 rtol = 0
end
