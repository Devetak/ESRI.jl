using Test
using ESRI
using SparseArrays
using Random

@testset "New Options: Combination Modes" begin
    Random.seed!(1001)
    N = 10
    W = sprand(N, N, 0.2) + 0.1I
    industry_ids = rand(1:3, N)
    essential_industry = [true, false, true]
    info = IndustryInfo(industry_ids, essential_industry)
    econ = ESRIEconomy(W, info)

    # Base run
    res_default = esri(econ, 1; details=true)
    
    # Check :min (should be default)
    res_min = esri(econ, 1; combine=:min, details=true)
    @test res_min.esri ≈ res_default.esri

    # Check :upstream
    res_up = esri(econ, 1; combine=:upstream, details=true)
    # Manual calculation from details
    expected_up = sum(econ.row_sums .* (1 .- res_up.upstream)) / econ.total_output
    @test res_up.esri ≈ expected_up

    # Check :downstream
    res_down = esri(econ, 1; combine=:downstream, details=true)
    # Manual calculation from details
    expected_down = sum(econ.row_sums .* (1 .- res_down.downstream)) / econ.total_output
    @test res_down.esri ≈ expected_down
    
    # Economy-wide combine
    scores_up = esri(econ; combine=:upstream)
    scores_down = esri(econ; combine=:downstream)
    @test scores_up[1] ≈ res_up.esri
    @test scores_down[1] ≈ res_down.esri
end

@testset "New Options: Custom Weights" begin
    Random.seed!(1002)
    N = 10
    W = sprand(N, N, 0.2) + 0.1I
    industry_ids = rand(1:3, N)
    essential_industry = [true, false, true]
    info = IndustryInfo(industry_ids, essential_industry)
    econ = ESRIEconomy(W, info)

    custom_weights = rand(N)
    
    # Single firm
    res_custom = esri(econ, 2; final_weights=custom_weights, details=true)
    expected = sum(custom_weights .* (1 .- min.(res_custom.upstream, res_custom.downstream))) / econ.total_output
    @test res_custom.esri ≈ expected

    # Economy wide
    scores_custom = esri(econ; final_weights=custom_weights)
    @test scores_custom[2] ≈ res_custom.esri
    @test esri(econ; final_weights = Float32.(custom_weights))[2] ≈ res_custom.esri

    # Weight length validation
    @test_throws DimensionMismatch esri(econ; final_weights=rand(N-1))
end

@testset "New Options: Custom Shock" begin
    Random.seed!(1003)
    N = 15
    W = sprand(N, N, 0.1) + 0.1I
    industry_ids = rand(1:3, N)
    essential_industry = [true, false, true]
    info = IndustryInfo(industry_ids, essential_industry)
    econ = ESRIEconomy(W, info)

    # All ones shock (no shock)
    my_shock_ones = ones(N)
    res_ones = esri(econ, 1; shock=my_shock_ones, details=true)
    @test res_ones.esri ≈ 0.0 atol=1e-10
    @test all(res_ones.upstream .≈ 1.0)
    @test all(res_ones.downstream .≈ 1.0)

    # Custom partial shock
    my_shock = rand(N)
    res_custom = esri(econ, 1; shock=my_shock, details=true)
    @test all(res_custom.upstream .<= my_shock .+ 1e-12)
    @test all(res_custom.downstream .<= my_shock .+ 1e-12)

    # Validation
    @test_throws DimensionMismatch esri(econ, 1; shock=rand(N-1))
    @test_throws DomainError esri(econ, 1; shock=fill(1.5, N))
    @test_throws DomainError esri(econ, 1; shock=fill(-0.1, N))

    # Test dedicated esri_shock API
    res_shock = esri_shock(econ, my_shock; details=true)
    @test res_shock.esri ≈ res_custom.esri atol=1e-12 rtol=0
    @test res_shock.upstream ≈ res_custom.upstream atol=1e-12 rtol=0
    @test res_shock.downstream ≈ res_custom.downstream atol=1e-12 rtol=0

    # compute_esri_shock wrapper
    s3 = compute_esri_shock(W, info, my_shock)
    @test s3 ≈ res_shock.esri
end

@testset "New Options: compute_esri overloads" begin
    Random.seed!(1004)
    N = 10
    W = sprand(N, N, 0.2) + 0.1I
    industry_ids = rand(1:3, N)
    essential_industry = [true, false, true]
    info = IndustryInfo(industry_ids, essential_industry)
    
    custom_weights = rand(N)
    
    # Economy wide
    s1 = compute_esri(W, info; final_weights=custom_weights, combine=:upstream)
    
    # Single firm
    s2 = compute_esri(W, info, 3; final_weights=custom_weights, combine=:upstream)
    
    @test s1[3] ≈ s2
end

@testset "New Options: Consistency with threads" begin
    Random.seed!(1005)
    N = 20
    W = sprand(N, N, 0.1) + 0.1I
    industry_ids = rand(1:3, N)
    essential_industry = [true, false, true]
    info = IndustryInfo(industry_ids, essential_industry)
    econ = ESRIEconomy(W, info)
    
    custom_weights = rand(N)
    
    s_serial = esri(econ; final_weights=custom_weights, combine=:downstream, threads=false)
    s_thread = esri(econ; final_weights=custom_weights, combine=:downstream, threads=true)
    
    @test s_serial ≈ s_thread
end
