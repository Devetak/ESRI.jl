using Test
using ESRIcascade

include("../benchmark/sparse_powerlaw_esri.jl")

function main(args = ARGS)
    n = isempty(args) ? 5_000 : parse(Int, args[1])
    result = benchmark_once(n; threaded = false, maxiter = 20, tol = 1e-3)
    print_benchmark(result)
    @test 0.0 <= first(result.score_range) <= last(result.score_range) <= 1.0 + 1e-6
end

main()
