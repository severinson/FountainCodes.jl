using FountainCodes, Test, Random
println("Starting tests")
@time @testset "GF256" begin include("GF256_test.jl") end
@time @testset "Gray sequence" begin include("Gray_test.jl") end
@time @testset "QMatrix" begin include("QMatrix_test.jl") end
@time @testset "Decoder" begin include("Decoder_test.jl") end
@time @testset "R10" begin include("R10_test.jl") end
@time @testset "RQ" begin include("RQ_test.jl") end
@time @testset "LT" begin include("LT_test.jl") end