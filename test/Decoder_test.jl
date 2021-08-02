# Basic decoder tests

using FountainCodes, LinearAlgebra, SparseArrays, Test

function rand_nonzero(rng::AbstractRNG, T::Type)
    rv = rand(rng, T)
    while iszero(rv)
        rv = rand(rng, T)    
    end
    rv
end

function rand_nonzero(rng::AbstractRNG, T::Type, dims...)
    rv = zeros(T, dims...)
    for i in 1:length(rv)
        rv[i] = rand_nonzero(rng, T)
    end
    rv
end

function test_decoder(A, b; expected_inactivations=nothing)
    x = transpose(A)*b
    println("Input source: $(Int.(b))")
    println("Input symbols: $(Int.(x))")
    tc = TestConstraints(A, x)
    decoder = Decoder(A)
    dec = decode(A, tc; decoder)
    # @show Int.(tc.A')
    println("Input source: $(Int.(b))")
    println("Input symbols: $(Int.(x))")
    println("Decoded symbols: $(Int.(dec))")    
    @test dec == b
    if !isnothing(expected_inactivations)
        @test decoder.metrics["inactivations"] == expected_inactivations
    end
    true
end

"""Test decoding with a diagonal constraint matrix."""
function test_diagonal(K::Integer, Tv=GF256)
    rng = MersenneTwister(123)
    constraints = [sparsevec([i], [rand_nonzero(rng, Tv)], K) for i in 1:K]
    A = hcat(constraints...)
    b = rand(rng, Tv, K)
    test_decoder(A, b, expected_inactivations=0)
end
# for K in [1, 10, 100, 200, 250, 254]
#     @test test_diagonal(Bool, K)
# end
for K in [1, 10, 100, 200, 250, 254]
    @test test_diagonal(K)
end

function test_bidiagonal(K::Integer, Tv=GF256)
    rng = MersenneTwister(123)
    dv = rand_nonzero(rng, Tv, K)
    ev = rand_nonzero(rng, Tv, K-1)
    A = sparse(Bidiagonal(dv, ev, :U))
    b = rand(rng, Tv, K)
    test_decoder(A, b, expected_inactivations=0)
end
for K in [10, 100, 200, 250, 254]
    @test test_bidiagonal(K)
end

function test_bidiagonal_permuted(K::Integer, Tv=GF256)
    rng = MersenneTwister(123)
    dv = rand_nonzero(rng, Tv, K)
    ev = rand_nonzero(rng, Tv, K-1)
    p = randperm(rng, K)
    q = randperm(rng, K)
    A = sparse(Bidiagonal(dv, ev, :U)[p, q])
    b = rand(rng, Tv, K)
    test_decoder(A, b, expected_inactivations=0)
end
for K in [10, 100, 200, 250, 254]
    @test test_bidiagonal(K)
end

function test_tridiagonal(K::Integer, Tv=GF256)
    rng = MersenneTwister(123)

    du = rand_nonzero(rng, Tv, K-1)
    d = rand_nonzero(rng, Tv, K)
    dl = rand_nonzero(rng, Tv, K-1)
    A = sparse(Tridiagonal(dl, d, du))
    b = rand(rng, Tv, K)
    test_decoder(A, b, expected_inactivations=1)
end
for K in [10, 100, 200, 250, 254]
    @test test_tridiagonal(K)
end

"""Test decoding for a dense constraint matrix."""
function test_dense(m::Integer, n::Integer=m, Tv=GF256)
    rng = MersenneTwister(123)    
    A = sparse(rand_nonzero(rng, Tv, m, n))
    b = rand(rng, Tv, m)
    test_decoder(A, b, expected_inactivations=m-1)
end
for K in [10, 100, 200, 250, 254]
    @test test_dense(K)
end
