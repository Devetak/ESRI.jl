"""
    compute_esri(weight_matrix, info::IndustryInfo; maxiter=100, tol=1e-2, verbose=false, threads=false)

Compute ESRI (Economic Systemic Risk Index) for each firm.

If `threads=true`, the outer firm loop is parallelized with `Threads.@threads`. Each worker allocates its own work buffers to avoid races; dense linear algebra inside may still use BLAS threads depending on your Julia/BLAS configuration.
"""
function compute_esri(
    weight_matrix,
    info::IndustryInfo;
    maxiter::Int = 100,
    tol::Real = 1e-2,
    verbose::Bool = false, threads::Bool = false, firm_indices::Union{Nothing,AbstractVector{<:Integer}} = nothing,
)
    n = length(info)  # number of firms
    esri = zeros(eltype(weight_matrix), n)

    # Precompute impact matrices
    upstream_impact_matrix = create_upstream_impact_matrix(weight_matrix)
    # Use CSR for sparse matrices, regular matrix for dense
    if weight_matrix isa SparseArrays.SparseMatrixCSC
        downstream_impact_matrix = compute_downstream_impact_matrix_csr(weight_matrix, info)
    else
        downstream_impact_matrix = compute_downstream_impact_matrix(weight_matrix, info)
    end

    # Precompute column and row sums
    column_sums = vec(sum(weight_matrix, dims = 1))
    row_sums = vec(sum(weight_matrix, dims = 2))
    total_output = sum(column_sums)

    firm_sel = firm_indices === nothing ? (1:n) : collect(firm_indices)
    if firm_indices !== nothing
        @assert all(1 .<= firm_sel .<= n) "firm_indices must be in 1:n"
        @assert length(unique(firm_sel)) == length(firm_sel) "firm_indices must be unique"
    end

    use_threads = threads && Threads.nthreads() > 1
    if use_threads
        if verbose
            @warn "Ignoring `verbose=true` because progress UI is disabled in threaded mode."
        end
        nt = Threads.nthreads()
        tid0 = Threads.threadid()
        Tlocal = eltype(weight_matrix)
        oneT = one(Tlocal)
        zeroT = zero(Tlocal)
        num_inds = num_industries(info)

        current_upstream = [Vector{Tlocal}(undef, n) for _ in 1:nt]
        previous_upstream = [Vector{Tlocal}(undef, n) for _ in 1:nt]
        visited = [Vector{Bool}(undef, n) for _ in 1:nt]
        current_downstream = [Vector{Tlocal}(undef, n) for _ in 1:nt]
        previous_downstream = [Vector{Tlocal}(undef, n) for _ in 1:nt]
        sigmas = [Vector{Tlocal}(undef, n) for _ in 1:nt]
        product_matrix = [Matrix{Tlocal}(undef, n, num_inds) for _ in 1:nt]
        temp_sums = [Vector{Tlocal}(undef, num_inds) for _ in 1:nt]
        final_shock = [Vector{Tlocal}(undef, n) for _ in 1:nt]

        Threads.@threads for k in eachindex(firm_sel)
            firm_idx = firm_sel[k]
            tid = Threads.threadid() - tid0
            cu = current_upstream[tid]
            pu = previous_upstream[tid]
            v = visited[tid]
            cd = current_downstream[tid]
            pd = previous_downstream[tid]
            sg = sigmas[tid]
            pm = product_matrix[tid]
            ts = temp_sums[tid]
            fs = final_shock[tid]

            # Upstream shock
            if column_sums[firm_idx] == zeroT
                fill!(cu, oneT)
                cu[firm_idx] = zeroT
            else
                upstream_shock!(
                    upstream_impact_matrix,
                    firm_idx,
                    cu,
                    pu,
                    v;
                    maxiter = maxiter,
                    tol = tol,
                    verbose = false,
                )
            end

            # Downstream shock
            if row_sums[firm_idx] == zeroT
                fill!(cd, oneT)
                cd[firm_idx] = zeroT
            else
                downstream_shock!(
                    downstream_impact_matrix,
                    info,
                    column_sums,
                    firm_idx,
                    sg,
                    pm,
                    cd,
                    pd,
                    ts;
                    maxiter = maxiter,
                    tol = tol,
                    verbose = false,
                )
            end

            # ESRI contribution for this firm
            acc = zeroT
            @inbounds for j in 1:n
                f = min(cu[j], cd[j])
                fs[j] = f
                acc += row_sums[j] * (oneT - f)
            end
            esri[firm_idx] = acc
        end
    end

    # Preallocate working vectors
    if !use_threads
        current_upstream = Vector{eltype(weight_matrix)}(undef, n)
        previous_upstream = Vector{eltype(weight_matrix)}(undef, n)
        visited = Vector{Bool}(undef, n)
        current_downstream = Vector{eltype(weight_matrix)}(undef, n)
        previous_downstream = Vector{eltype(weight_matrix)}(undef, n)
        sigmas = Vector{eltype(weight_matrix)}(undef, n)
        product_matrix = Matrix{eltype(weight_matrix)}(undef, n, num_industries(info))
        temp_sums = Vector{eltype(weight_matrix)}(undef, num_industries(info))
        final_shock = Vector{eltype(weight_matrix)}(undef, n)

        # Progress bar if verbose
        iter_range = verbose ? ProgressBar(firm_sel; total = length(firm_sel)) : firm_sel

        # Compute ESRI for each firm
        for firm_idx in iter_range
            # Handle upstream shock
            if column_sums[firm_idx] == 0.0
                fill!(current_upstream, 1.0)
                current_upstream[firm_idx] = 0.0
            else
                upstream_shock!(
                    upstream_impact_matrix,
                    firm_idx,
                    current_upstream,
                    previous_upstream,
                    visited;
                    maxiter = maxiter,
                    tol = tol,
                    verbose = false,
                )
            end

            # Handle downstream shock
            if row_sums[firm_idx] == 0.0
                fill!(current_downstream, 1.0)
                current_downstream[firm_idx] = 0.0
            else
                downstream_shock!(
                    downstream_impact_matrix,
                    info,
                    column_sums,
                    firm_idx,
                    sigmas,
                    product_matrix,
                    current_downstream,
                    previous_downstream,
                    temp_sums;
                    maxiter = maxiter,
                    tol = tol,
                    verbose = false,
                )
            end

            # Compute final shock as minimum of upstream and downstream
            @inbounds for j = 1:n
                final_shock[j] = min(current_upstream[j], current_downstream[j])
            end

            # Compute ESRI contribution
            esri[firm_idx] = 0.0
            @inbounds for j = 1:n
                esri[firm_idx] += row_sums[j] * (1.0 - final_shock[j])
            end
        end
    end

    # Normalize by total output
    if total_output > 0
        esri ./= total_output
    end

    return esri
end
