# Basic decoder tests

using FountainCodes, LinearAlgebra, SparseArrays, Test

function test_decoder(A, b; expected_inactivations=nothing)
    x = transpose(A)*b
    # println("Input source: $(Int.(b))")
    # println("Input symbols: $(Int.(x))")
    tc = TestConstraints(A, x)
    decoder = Decoder(A)
    dec = decode(A, tc; decoder)
    # @show Int.(tc.A')
    # println("Input source: $(Int.(b))")
    # println("Input symbols: $(Int.(x))")
    # println("Decoded symbols: $(Int.(dec))")    
    if eltype(A) <: AbstractFloat
        @test dec ≈ b
    else
        @test dec == b
    end
    if !isnothing(expected_inactivations)
        @test decoder.metrics["inactivations"] == expected_inactivations
    end
    true
end

"""Test decoding with a diagonal constraint matrix."""
function test_diagonal(K::Integer; Tv=GF256)
    rng = MersenneTwister(123)
    constraints = [sparsevec([i], [rand_nonzero(rng, Tv)], K) for i in 1:K]
    A = hcat(constraints...)
    b = rand(rng, Tv, K)
    test_decoder(A, b, expected_inactivations=0)
end
for K in [1, 10, 100, 200, 250, 254]
    test_diagonal(K)
end

function test_bidiagonal(K::Integer; Tv=GF256, Ti=Int, permuted=false)
    rng = MersenneTwister(123)
    dv = rand_nonzero(rng, Tv, K)
    ev = rand_nonzero(rng, Tv, K-1)
    b = rand(rng, Tv, K)    
    if permuted
        p = randperm(rng, K)
        q = randperm(rng, K)        
        A = SparseMatrixCSC{Tv,Ti}(sparse(Bidiagonal(dv, ev, :U)[p, q]))
    else
        A = SparseMatrixCSC{Tv,Ti}(sparse(Bidiagonal(dv, ev, :U)))
    end
    test_decoder(A, b, expected_inactivations=0)
end
for K in [10, 100, 200, 250, 254]
    test_bidiagonal(K, permuted=false)
    test_bidiagonal(K, permuted=true)
end
test_bidiagonal(10, Ti=Int32, permuted=false)

function test_tridiagonal(K::Integer; Tv=GF256, permuted=false)
    rng = MersenneTwister(123)
    du = rand_nonzero(rng, Tv, K-1)
    d = rand_nonzero(rng, Tv, K)
    dl = rand_nonzero(rng, Tv, K-1)
    b = rand(rng, Tv, K)
    if permuted    
        p = randperm(rng, K)
        q = randperm(rng, K)                
        A = sparse(Tridiagonal(dl, d, du)[p, q])
    else
        A = sparse(Tridiagonal(dl, d, du))
    end
    test_decoder(A, b, expected_inactivations=1)
end
for K in [10, 100, 200, 250, 254]
    test_tridiagonal(K, permuted=false)
    test_tridiagonal(K, permuted=true)
end
test_tridiagonal(10, Tv=Float64, permuted=false)

"""Test decoding for a dense constraint matrix."""
function test_dense(m::Integer, n::Integer=m; Tv=GF256)
    rng = MersenneTwister(123)    
    A = sparse(rand_nonzero(rng, Tv, m, n))
    b = rand(rng, Tv, m)
    test_decoder(A, b, expected_inactivations=m-1)
end
for K in [10, 100, 200, 250, 254]
    test_dense(K)
end

# """Test decoding with each value consisting of an array"""
# function test_vector_values(m::Integer, n::Integer=m, dims=(1,); Tv=GF256)
#     rng = MersenneTwister(123)    
#     A = sparse(rand_nonzero(rng, Tv, m, n))
#     b = [rand(rng, Tv, dims...) for _ in 1:m]
#     test_decoder(A, b, expected_inactivations=m-1)
# end
# test_vector_values(10)

"""Test that connected components are selected correctly."""
function test_inactivations_1(Tv=GF256)
    rng = MersenneTwister(123)
    Is = [
        (1, 1), (1, 2),
        (2, 2), (2, 3),
        (3, 3), (3, 4),
        (4, 4), (4, 5),
        (5, 1), (5, 2),
        (6, 2), (6, 3),
        (7, 4), (7, 5)
    ]
    M = zeros(Tv, 7, 5)
    for (i, j) in Is
        M[i, j] = rand_nonzero(rng, Tv)
    end
    b = rand(rng, Tv, 5)    
    p = randperm(rng, size(M, 1))
    q = randperm(rng, size(M, 2))
    A = sparse(M[p, q]')
    test_decoder(A, b, expected_inactivations=1)
end
test_inactivations_1()

"""Test that connected components are selected correctly."""
function test_inactivations_2(Tv=GF256)
    rng = MersenneTwister(123)
    Is = [
        # Component 1 (4 symbols)
        (1, 1), (1, 2),
        (2, 2), (2, 3),
        (3, 3), (3, 4),
        # Component 2 (3 symbols)
        (4, 5), (4, 6),
        (5, 6), (5, 7),
    ]
    M = zeros(Tv, 8, 7)
    for (i, j) in Is
        M[i, j] = rand_nonzero(rng, Tv)
    end
    M[6, 1:4] .= rand_nonzero(rng, Tv, 4)
    M[7, 1:4] .= rand_nonzero(rng, Tv, 4)
    M[8, 1:4] .= rand_nonzero(rng, Tv, 4)
    M[6, 5] = rand_nonzero(rng, Tv)    
    M[7, 6] = rand_nonzero(rng, Tv)    
    M[8, 7] = rand_nonzero(rng, Tv)    
    b = rand(rng, Tv, 7)
    p = randperm(rng, size(M, 1))
    q = randperm(rng, size(M, 2))
    A = sparse(M[p, q]')
    test_decoder(A, b, expected_inactivations=1)
end
test_inactivations_2()