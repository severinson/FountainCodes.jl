
@testset "Properties" begin
    @test iszero(GF256(0))
    @test !iszero(GF256(1))
    @test !iszero(GF256(255))
    @test isone(GF256(1))
    @test !isone(GF256(0))
    @test !isone(GF256(255))
    @test zero(GF256) == GF256(0)
    @test zero(GF256) == zero(GF256(1))
    @test one(GF256) == GF256(1)
    @test one(GF256) == one(GF256(0))    
end

@testset "Logarithms and exponentiation" begin
    a = GF256(0)
    @test log(exp(a)) == a
    for x in 1:254
        a = GF256(x)
        @test log(exp(a)) == a
        @test exp(log(a)) == a                
    end
    a = GF256(255)
    @test exp(log(a)) == a
end

@testset "Division and multiplication" begin
    β = GF256(43) # chosen arbitrarily
    @test β * zero(GF256) == zero(GF256)
    @test β * one(GF256) == β
    for x in 1:255
        @test (β * GF256(x)) / GF256(x) == β
        @test (β / GF256(x)) * GF256(x) == β
    end
end

@testset "Addition and subtraction" begin
    β = GF256(43) # chosen arbitrarily 
    @test β + zero(GF256) == β
    @test β - zero(GF256) == β
    for x in 0:255
        a = GF256(x)
        @test iszero(a+a)
        @test iszero(a-a)
    end
end

# TODO: remove
# @test GF256(2) + GF256(10) - GF256(10) == GF256(2)
# @test GF256(2) * GF256(10) / GF256(10) == GF256(2)
# @test GF256(1) + GF256(2) * GF256(10) - GF256(10) - GF256(10) == GF256(21)
# @test GF256(13) * GF256(1) == GF256(13)
# @test GF256(13) / GF256(1) == GF256(13)
# @test GF256(0) * GF256(1) == GF256(0)
# @test GF256(0) / GF256(1) == GF256(0)
# @test [GF256(i) for i in 0:255].*zero(GF256) == zeros(GF256, 256)
# @test [GF256(i) for i in 0:255].*one(GF256) == [GF256(i) for i in 0:255]
# @test [GF256(i) for i in 0:255].*GF256(2) == [GF256(2)*GF256(i) for i in 0:255]
# # @test GF256(11) + false == GF256(11)
# # @test GF256(11) + true == GF256(11) + GF256(1)
# # @test GF256(11) * false == GF256(0)
# # @test GF256(11) * true == GF256(11)
# # @test GF256(11) / true == GF256(11)
# @test GF256(123)^0 == GF256(1)
# @test GF256(123)^1 == GF256(123)
# @test GF256(123)^2 == GF256(123) * GF256(123)
# @test GF256(46)^3 == GF256(46) * GF256(46) * GF256(46)
# @test GF256(10)^20 == reduce(*, (GF256(10) for _ in 1:20))
# @test GF256(56)^30 == reduce(*, (GF256(56) for _ in 1:30))