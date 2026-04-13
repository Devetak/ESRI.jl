using Test
using Random
using SparseArrays
using LinearAlgebra: I, BLAS

using ESRI

include("test_helpers.jl")
include("test_types.jl")
include("test_primitives.jl")
include("test_core_api.jl")
include("test_new_options.jl")
include("test_regressions.jl")

if "quality" in ARGS
    include("test_quality.jl")
end
