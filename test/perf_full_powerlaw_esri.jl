using ESRIcascade
import ESRIcascade as ESRI

include("../benchmark/sparse_powerlaw_esri.jl")

function main(args = ARGS)
    n = isempty(args) ? 5_000 : parse(Int, args[1])
    modes = length(args) < 2 ? (:truncated_tail, :heavy_tail) : (Symbol(args[2]),)

    for mode in modes
        result = benchmark_once(n; mode = mode, threaded = false, maxiter = 20, tol = 1e-3)
        print_benchmark(result)
    end
end

main()
