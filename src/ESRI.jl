module ESRI

using LinearAlgebra
using SparseArrays
using SparseMatricesCSR
using ProgressBars

include("types.jl")
include("impact.jl")
include("upstream.jl")
include("downstream.jl")
include("esri_computation.jl")

export IndustryInfo, ESRIEconomy, ESRIResult, esri, compute_esri

end
