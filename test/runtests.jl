#!/usr/bin/env julia

using RaptorCodes, Base.Test
println("Starting tests")
# @time @testset "Gaussian elimination" begin include("solve_test.jl") end
@time @testset "Gray sequence" begin include("gray_test.jl") end
@time @testset "Raptor encoder" begin include("R10Encode_test.jl") end
@time @testset "Raptor decoder" begin include("decoder_test.jl") end
