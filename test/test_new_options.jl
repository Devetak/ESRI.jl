using Test
using ESRI
using SparseArrays

function deterministic_fixture()
    # 3-firm hand-solvable network:
    # firm 1 supplies firms 1 and 2; firm 2 self-supplies; firm 3 isolated self-supply.
    W = [
        1.0 1.0 0.0
        0.0 1.0 0.0
        0.0 0.0 1.0
    ]
    info = IndustryInfo([1, 1, 1], [true])
    econ = ESRIEconomy(W, info)
    return W, info, econ
end

@testset "New Options: Combination Modes (hand-computed)" begin
    _, _, econ = deterministic_fixture()
    maxiter = 100

    # Shock firm 1 only. For this network, upstream reaches [0,1,1] exactly.
    # Downstream iteration for firm 2 follows x_{k+1} = x_k / (1 + x_k), x_0 = 1.
    # Under this solver's update order, after `maxiter = m` iterations the observed value is x = 1 / m.
    details = esri(econ, 1; details = true, maxiter = maxiter, tol = 1e-12)
    x = 1 / maxiter
    @test details.upstream ≈ [0.0, 1.0, 1.0] atol = 1e-12 rtol = 0
    @test details.downstream ≈ [0.0, x, 1.0] atol = 1e-12 rtol = 0

    # row_sums = [2, 1, 1], total_output = 4
    # :min => (2*(1-0) + 1*(1-x) + 1*(1-1))/4 = (3-x)/4
    # :upstream => (2*(1-0) + 1*(1-1) + 1*(1-1))/4 = 1/2
    # :downstream => same as :min here
    expected_min = (3 - x) / 4
    @test esri(econ, 1; combine = :min, maxiter = maxiter, tol = 1e-12) ≈ expected_min atol = 1e-12 rtol = 0
    @test esri(econ, 1; combine = :upstream, maxiter = maxiter, tol = 1e-12) ≈ (1 / 2) atol = 1e-12 rtol = 0
    @test esri(econ, 1; combine = :downstream, maxiter = maxiter, tol = 1e-12) ≈ expected_min atol = 1e-12 rtol = 0
    @test details.esri ≈ expected_min atol = 1e-12 rtol = 0

    scores_up = esri(econ; combine = :upstream, maxiter = maxiter, tol = 1e-12)
    scores_down = esri(econ; combine = :downstream, maxiter = maxiter, tol = 1e-12)
    @test scores_up[1] ≈ (1 / 2) atol = 1e-12 rtol = 0
    @test scores_down[1] ≈ expected_min atol = 1e-12 rtol = 0
end

@testset "New Options: Custom Weights (hand-computed)" begin
    _, _, econ = deterministic_fixture()
    custom_weights = [5.0, 7.0, 11.0]
    maxiter = 100

    # After `maxiter = m`, min(upstream, downstream) is [0, 1/m, 1].
    # weighted loss = 5*(1-0) + 7*(1-1/m) + 11*(1-1)
    x = 1 / maxiter
    expected = (5 + 7 * (1 - x)) / 4
    res_custom = esri(econ, 1; final_weights = custom_weights, details = true, maxiter = maxiter, tol = 1e-12)
    @test res_custom.esri ≈ expected atol = 1e-12 rtol = 0

    scores_custom = esri(econ; final_weights = custom_weights, maxiter = maxiter, tol = 1e-12)
    @test scores_custom[1] ≈ expected atol = 1e-12 rtol = 0
    @test esri(econ; final_weights = Float32.(custom_weights), maxiter = maxiter, tol = 1e-12)[1] ≈ expected atol = 1e-6 rtol = 0

    @test_throws DimensionMismatch esri(econ; final_weights = [1.0, 2.0])
end

@testset "New Options: Custom Shock (deterministic)" begin
    W, info, econ = deterministic_fixture()

    shock_ones = ones(3)
    res_ones = esri(econ, 1; shock = shock_ones, details = true, maxiter = 100, tol = 1e-12)
    @test res_ones.esri ≈ 0.0 atol = 1e-12
    @test res_ones.upstream == ones(3)
    @test res_ones.downstream == ones(3)

    shock_zeros = zeros(3)
    res_zeros = esri_shock(econ, shock_zeros; details = true, maxiter = 100, tol = 1e-12)
    @test res_zeros.esri ≈ 1.0 atol = 1e-10
    @test all(isapprox.(res_zeros.upstream, 0.0; atol = 1e-12, rtol = 0))
    @test all(isapprox.(res_zeros.downstream, 0.0; atol = 1e-12, rtol = 0))

    shock_partial = [0.2, 0.8, 1.0]
    res_single = esri(econ, 1; shock = shock_partial, details = true, maxiter = 100, tol = 1e-12)
    @test all(res_single.upstream .<= shock_partial .+ 1e-12)
    @test all(res_single.downstream .<= shock_partial .+ 1e-12)

    @test_throws DimensionMismatch esri(econ, 1; shock = [0.1, 0.9])
    @test_throws DomainError esri(econ, 1; shock = [1.2, 0.8, 0.5])
    @test_throws DomainError esri(econ, 1; shock = [-0.1, 0.8, 0.5])

    res_shock = esri_shock(econ, shock_partial; details = true, maxiter = 100, tol = 1e-12)
    @test res_shock.esri ≈ res_single.esri atol = 1e-12 rtol = 0
    @test res_shock.upstream ≈ res_single.upstream atol = 1e-12 rtol = 0
    @test res_shock.downstream ≈ res_single.downstream atol = 1e-12 rtol = 0

    wrapped = compute_esri_shock(W, info, shock_partial; maxiter = 100, tol = 1e-12)
    @test wrapped ≈ res_shock.esri atol = 1e-12 rtol = 0

    # components branches
    res_up_comp = esri_shock(econ, shock_partial; components = :upstream, maxiter = 100, tol = 1e-12)
    @test haskey(res_up_comp, :esri)
    @test haskey(res_up_comp, :upstream)
    @test !haskey(res_up_comp, :downstream)
    @test res_up_comp.upstream ≈ res_shock.upstream atol = 1e-12 rtol = 0
    @test res_up_comp.esri ≈ res_shock.esri atol = 1e-12 rtol = 0

    res_down_comp = esri_shock(econ, shock_partial; components = :downstream, maxiter = 100, tol = 1e-12)
    @test haskey(res_down_comp, :esri)
    @test haskey(res_down_comp, :downstream)
    @test !haskey(res_down_comp, :upstream)
    @test res_down_comp.downstream ≈ res_shock.downstream atol = 1e-12 rtol = 0
    @test res_down_comp.esri ≈ res_shock.esri atol = 1e-12 rtol = 0

    # Explicit combine formula checks
    res_up_only = esri_shock(econ, shock_partial; components = :upstream, combine = :upstream, maxiter = 100, tol = 1e-12)
    expected_up = sum(econ.row_sums .* (1 .- res_up_only.upstream)) / econ.total_output
    @test res_up_only.esri ≈ expected_up atol = 1e-12 rtol = 0

    res_down_only = esri_shock(econ, shock_partial; components = :downstream, combine = :downstream, maxiter = 100, tol = 1e-12)
    expected_down = sum(econ.row_sums .* (1 .- res_down_only.downstream)) / econ.total_output
    @test res_down_only.esri ≈ expected_down atol = 1e-12 rtol = 0
end

@testset "New Options: compute_esri overloads (deterministic)" begin
    W, info, _ = deterministic_fixture()
    custom_weights = [3.0, 2.0, 1.0]

    economy_scores = compute_esri(W, info; final_weights = custom_weights, combine = :upstream, maxiter = 100, tol = 1e-12)
    single_score = compute_esri(W, info, 1; final_weights = custom_weights, combine = :upstream, maxiter = 100, tol = 1e-12)

    @test economy_scores[1] ≈ single_score atol = 1e-12 rtol = 0
end

@testset "New Options: Consistency with threads (deterministic)" begin
    W = sparse(
        [1, 1, 2, 2, 3, 4, 5, 6, 3, 4],
        [1, 2, 2, 3, 3, 4, 5, 6, 4, 5],
        [1.0, 0.4, 1.0, 0.2, 1.0, 1.0, 1.0, 1.0, 0.3, 0.1],
        6,
        6,
    )
    info = IndustryInfo([1, 1, 2, 2, 1, 2], [true, false])
    econ = ESRIEconomy(W, info)
    custom_weights = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]

    s_serial = esri(econ; final_weights = custom_weights, combine = :downstream, threads = false, maxiter = 80, tol = 1e-10)
    s_thread = esri(econ; final_weights = custom_weights, combine = :downstream, threads = true, maxiter = 80, tol = 1e-10)
    @test s_serial ≈ s_thread atol = 1e-12 rtol = 0
end
