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
