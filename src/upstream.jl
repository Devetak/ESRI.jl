function upstream_step!(
    next_step::AbstractVector,
    impact_matrix::AbstractMatrix,
    previous_step::AbstractVector,
)
    mul!(next_step, transpose(impact_matrix), previous_step)
    return next_step
end

function upstream_step!(
    next_step::AbstractVector,
    impact_matrix::AbstractMatrix,
    previous_step::AbstractVector,
    visited::AbstractVector{Bool},
)
    # For dense matrices, visited is not needed, just call the plain version
    upstream_step!(next_step, impact_matrix, previous_step)
    return next_step
end

function upstream_step!(
    next_step::AbstractVector,
    impact_matrix::SparseMatrixCSC,
    previous_step::AbstractVector,
    visited::AbstractVector{Bool},  # kept for API, unused
)
    rows = impact_matrix.rowval
    vals = impact_matrix.nzval
    colptr = impact_matrix.colptr
    ncols = size(impact_matrix, 2)

    zeroT = zero(eltype(next_step))
    oneT = one(eltype(next_step))

    @inbounds for col = 1:ncols
        start_idx = colptr[col]
        stop_idx = colptr[col+1] - 1

        if start_idx > stop_idx
            # structurally empty column
            next_step[col] = oneT
            continue
        end

        acc = zeroT
        @simd for idx = start_idx:stop_idx
            row = rows[idx]
            acc = muladd(vals[idx], previous_step[row], acc)
        end
        next_step[col] = acc
    end

    return next_step
end

@inline function _linf_distance(a::AbstractVector, b::AbstractVector)
    T = promote_type(eltype(a), eltype(b), Float64)
    max_diff = zero(T)
    @inbounds @simd for i in eachindex(a, b)
        diff = abs(a[i] - b[i])
        if diff > max_diff
            max_diff = diff
        end
    end
    return max_diff
end

function upstream_shock!(
    impact_matrix::AbstractMatrix,
    firm_idx::Integer,
    current_upstream::AbstractVector,
    previous_upstream::AbstractVector,
    visited::AbstractVector{Bool};
    maxiter::Int = 100,
    tol::Real = 1e-2,
    verbose::Bool = false,
    every::Int = 10,
)
    oneT = one(eltype(current_upstream))
    zeroT = zero(eltype(current_upstream))

    current_upstream .= oneT
    previous_upstream .= oneT
    current_upstream[firm_idx] = zeroT
    previous_upstream[firm_idx] = zeroT

    # local aliases that we can swap
    x = current_upstream
    y = previous_upstream

    for iter = 1:maxiter
        upstream_step!(x, impact_matrix, y, visited)
        x[firm_idx] = zeroT

        distance = _linf_distance(x, y)
        if verbose && ((iter - 1) % every == 0)
            @info "upstream iteration" iter distance
        end
        if distance < tol
            # ensure the final result ends up in current_upstream
            if x !== current_upstream
                copyto!(current_upstream, x)
            end
            return nothing
        end

        # swap roles for next iteration
        x, y = y, x
    end

    # maxiter reached, y holds last iterate
    if y !== current_upstream
        copyto!(current_upstream, y)
    end
    return nothing
end

function upstream_shock!(
    impact_matrix::AbstractMatrix,
    firm_idx::Integer,
    current_upstream::AbstractVector,
    previous_upstream::AbstractVector;
    maxiter::Int = 100,
    tol::Real = 1e-2,
    verbose::Bool = false,
    every::Int = 10,
)
    visited = fill(false, length(current_upstream))
    upstream_shock!(
        impact_matrix,
        firm_idx,
        current_upstream,
        previous_upstream,
        visited;
        maxiter = maxiter,
        tol = tol,
        verbose = verbose,
        every = every,
    )
end
