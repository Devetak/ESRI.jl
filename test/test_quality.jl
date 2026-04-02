using Test
using Aqua
using ESRI

@testset "Aqua quality checks" begin
    Aqua.test_all(ESRI; ambiguities = false)
end
