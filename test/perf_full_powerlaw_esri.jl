using Test
using ESRI

include("../benchmark/sparse_powerlaw_esri.jl")

function main()
    for mode in (:truncated_tail, :heavy_tail)
        result = benchmark_once(5_000; mode = mode, threaded = false, maxiter = 20, tol = 1e-3)
        print_benchmark(result)
        @test 0.0 <= first(result.score_range) <= last(result.score_range) <= 1.0 + 1e-6
    end
end

main()
