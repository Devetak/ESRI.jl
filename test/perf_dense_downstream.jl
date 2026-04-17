using ESRIcascade
using Random
using Statistics

function old_accumulate_downstream_components!(
    essential_matrix::AbstractMatrix,
    nonessential_vector::AbstractVector,
    downstream_vector::AbstractVector,
    sigmas::AbstractVector,
    essential_impact::AbstractMatrix,
    nonessential_impact::AbstractMatrix,
    info::IndustryInfo,
    industry_map::AbstractVector{Int} = info.industry_of_firm,
)
    T = eltype(essential_matrix)
    oneT = one(T)

    fill!(essential_matrix, zero(T))
    fill!(nonessential_vector, zero(T))

    nfirms = length(downstream_vector)
    @inbounds for supplier = 1:nfirms
        industry = industry_map[supplier]
        factor = sigmas[supplier] * (oneT - downstream_vector[supplier])
        if factor == 0
            continue
        end
        for customer = 1:size(essential_matrix, 1)
            val = essential_impact[supplier, customer]
            if val != 0
                essential_matrix[customer, industry] += factor * val
            end
            val_nonessential = nonessential_impact[supplier, customer]
            if val_nonessential != 0
                nonessential_vector[customer] += factor * val_nonessential
            end
        end
    end
    return nothing
end

function median_elapsed(f; samples::Int = 9, reps::Int = 15)
    times = Vector{Float64}(undef, samples)
    for s in 1:samples
        t0 = time_ns()
        for _ in 1:reps
            f()
        end
        times[s] = (time_ns() - t0) / 1e9
    end
    return median(times)
end

function main()
    Random.seed!(44)
    n = 280
    num_industries = 12

    W = rand(n, n)
    info = IndustryInfo(rand(1:num_industries, n), rand(Bool, num_industries))
    essential_impact, nonessential_impact = ESRIcascade.compute_downstream_impact_matrices(W, info)

    downstream = rand(n)
    sigmas_base = rand(n)

    em_old = zeros(n, num_industries)
    nv_old = zeros(n)
    em_new = zeros(n, num_industries)
    nv_new = zeros(n)
    sigmas = similar(sigmas_base)

    old_call! = () -> begin
        copyto!(sigmas, sigmas_base)
        old_accumulate_downstream_components!(
            em_old,
            nv_old,
            downstream,
            sigmas,
            essential_impact,
            nonessential_impact,
            info,
        )
    end

    new_call! = () -> begin
        copyto!(sigmas, sigmas_base)
        ESRIcascade._accumulate_downstream_components!(
            em_new,
            nv_new,
            downstream,
            sigmas,
            essential_impact,
            nonessential_impact,
            info,
        )
    end

    old_call!()
    new_call!()

    if !(isapprox(em_old, em_new; atol = 1e-12, rtol = 1e-12) &&
         isapprox(nv_old, nv_new; atol = 1e-12, rtol = 1e-12))
        error("dense accumulation mismatch between old and new implementation")
    end

    old_median = median_elapsed(old_call!)
    new_median = median_elapsed(new_call!)
    speedup = old_median / new_median

    println("old_median_s=", old_median)
    println("new_median_s=", new_median)
    println("speedup_x=", speedup)

    if new_median > old_median
        error("performance regression: new dense accumulation is slower")
    end
end

main()
