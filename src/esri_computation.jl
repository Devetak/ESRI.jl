function compute_esri(
    weight_matrix,
    info::IndustryInfo;
    maxiter::Int = 100,
    tol::Real = 1e-2,
    verbose::Bool = false,
)
    n = length(info)
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

    # Preallocate working vectors
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
    iter_range = verbose ? ProgressBar(1:n; total = n) : (1:n)

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
        @inbounds for j = 1:n
            esri[firm_idx] += row_sums[j] * (1.0 - final_shock[j])
        end
    end

    # Normalize by total output
    if total_output > 0
        esri ./= total_output
    end

    return esri
end
