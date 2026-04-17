using Test
using Aqua
using ESRIcascade

@testset "Aqua quality checks" begin
    Aqua.test_all(ESRIcascade; ambiguities = false)
end
