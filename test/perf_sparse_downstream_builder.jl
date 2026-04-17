using ESRIcascade
using Random
using SparseArrays
using SparseMatricesCSR
using Statistics

function old_sparse_downstream_builder(
    weight_matrix::SparseMatrixCSC{T},
    info::IndustryInfo,
) where {T<:Real}
    nrows, ncols = size(weight_matrix)
    rows = weight_matrix.rowval
    colptr = weight_matrix.colptr
    vals = weight_matrix.nzval
    essential_vals = zeros(T, length(vals))
    nonessential_vals = zeros(T, length(vals))

    num_inds = ESRIcascade.num_industries(info)
    essential_by_industry = zeros(T, num_inds)

    @inbounds for col = 1:ncols
        start_idx = colptr[col]
        stop_idx = colptr[col + 1] - 1
        if start_idx > stop_idx
            continue
        end

        fill!(essential_by_industry, zero(T))
        all_suppliers_total = zero(T)
        for idx = start_idx:stop_idx
            row = rows[idx]
            val = vals[idx]
            all_suppliers_total += val
            if ESRIcascade.is_essential(info, row)
                essential_by_industry[ESRIcascade.get_industry(info, row)] += val
            end
        end

        for idx = start_idx:stop_idx
            row = rows[idx]
            val = vals[idx]
            if ESRIcascade.is_essential(info, row)
                denom = essential_by_industry[ESRIcascade.get_industry(info, row)]
                essential_vals[idx] = denom == 0 ? zero(T) : val / denom
            else
                nonessential_vals[idx] = all_suppliers_total == 0 ? zero(T) : val / all_suppliers_total
            end
        end
    end

    essential_csc = SparseMatrixCSC(nrows, ncols, copy(colptr), copy(rows), essential_vals)
    nonessential_csc = SparseMatrixCSC(nrows, ncols, copy(colptr), copy(rows), nonessential_vals)

    I1, J1, V1 = findnz(essential_csc)
    I2, J2, V2 = findnz(nonessential_csc)
    return sparsecsr(I1, J1, V1, nrows, ncols), sparsecsr(I2, J2, V2, nrows, ncols)
end

function median_elapsed(f; samples::Int = 9, reps::Int = 10)
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
    Random.seed!(7)
    n = 1600
    density = 0.006
    inds = 25
    W = sprand(n, n, density)
    info = IndustryInfo(rand(1:inds, n), rand(Bool, inds))

    old_call = () -> old_sparse_downstream_builder(W, info)
    new_call = () -> ESRIcascade.compute_downstream_impact_matrices(W, info)

    old_ess, old_non = old_call()
    new_ess, new_non = new_call()
    if !(isapprox(old_ess.nzval, new_ess.nzval; atol = 1e-12, rtol = 1e-12) &&
         isapprox(old_non.nzval, new_non.nzval; atol = 1e-12, rtol = 1e-12))
        error("sparse downstream builder mismatch between old and new implementation")
    end

    old_median = median_elapsed(old_call)
    new_median = median_elapsed(new_call)
    speedup = old_median / new_median

    println("old_sparse_builder_median_s=", old_median)
    println("new_sparse_builder_median_s=", new_median)
    println("speedup_x=", speedup)

    if new_median > old_median
        error("sparse downstream builder regression: new path is slower")
    end
end

main()
