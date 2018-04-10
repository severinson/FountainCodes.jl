#!/usr/bin/env julia

using RaptorCodes, Base.Test
println("Starting tests")
@time @testset "Bounds tests" begin include("Bound_test.jl") end
@time @testset "Numerical function inversion" begin include("Numinv_test.jl") end
@time @testset "Finite field arithmetic" begin include("Arithmetic_test.jl") end
@time @testset "Soliton distribution" begin include("Soliton_test.jl") end
@time @testset "Gray sequence" begin include("gray_test.jl") end
@time @testset "R10 encoder" begin include("R10Encode_test.jl") end
@time @testset "Raptor decoder" begin include("decoder_test.jl") end
@time @testset "RQ encoder" begin include("RQ_test.jl") end
@time @testset "LT encoder" begin include("LT_test.jl") end
