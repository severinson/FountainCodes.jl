# Basic decoder tests

using FountainCodes, LinearAlgebra, SparseArrays, Test

"""Test decoding with a diagonal constraint matrix."""
function test_diagonal(CT, K::Integer)
    constraints = [sparsevec([i], [one(CT)], K) for i in 1:K]
    A = hcat(constraints...)
    Vs = [GF256(i%256) for i in 0:K-1]
    # d = Decoder{CT}(K)
    # dec = decode(constraints, Vs, decoder=d)
    decoder = Decoder(A)
    dec = decode(A, Vs; decoder)
    for i in 1:K
        if dec[i] != Vs[i]
            error("expected $(Vs[i]), but got $(dec[i])")
        end
    end
    if decoder.metrics["inactivations"] != 0
        error("expected 0 inactivations, but got $(decoder.metrics["inactivations"])")
    end
    true
end
# for K in [1, 10, 100, 200, 250, 254]
#     @test test_diagonal(Bool, K)
# end
for K in [1, 10, 100, 200, 250, 254]
    @test test_diagonal(GF256, K)
end

function test_bidiagonal(K::Integer, Tv=GF256)
    dv = ones(Tv, K)
    ev = ones(Tv, K-1)
    A = sparse(Bidiagonal(dv, ev, :U))
    Vs = [GF256(i%256) for i in 0:K-1]
    decoder = Decoder(A)
    dec = decode(A, Vs; decoder)
    for i in 1:K
        if dec[i] != Vs[i]
            error("expected $(Vs[i]), but got $(dec[i])")
        end
    end
    if decoder.metrics["inactivations"] != 0
        error("expected 0 inactivations, but got $(decoder.metrics["inactivations"])")
    end        
    true
end
for K in [10, 100, 200, 250, 254]
    @test test_bidiagonal(K)
end

function test_tridiagonal(K::Integer, Tv=GF256)
    du = ones(Tv, K-1)
    d = ones(Tv, K)
    dl = ones(Tv, K-1)
    A = sparse(Tridiagonal(dl, d, du))
    Vs = [GF256(i%256) for i in 0:K-1]
    decoder = Decoder(A)
    dec = decode(A, Vs; decoder)
    for i in 1:K
        if dec[i] != Vs[i]
            error("expected $(Vs[i]), but got $(dec[i])")
        end
    end
    if decoder.metrics["inactivations"] != 1
        error("expected 1 inactivation, but got $(decoder.metrics["inactivations"])")
    end        
    true
end
for K in [10, 100, 200, 250, 254]
    @test test_tridiagonal(K)
end

"""Test decoding for a dense constraint matrix."""
function test_dense(K::Integer)
    Is = collect(1:K)
    Cs = [GF256(i) for i in 1:K]
    constraints = [sparsevec(Is, circshift(Cs, i), K) for i in 0:K-1]
    A = hcat(constraints...)    
    src = [GF256(i%256) for i in 0:K-1]
    enc = [dot(constraint, src) for constraint in constraints]
    # d = Decoder{GF256}(K)
    # dec = decode(constraints, enc, decoder=d)
    decoder = Decoder(A)
    dec = decode(A, enc; decoder)
    for i in 1:K
        if dec[i] != src[i]
            error("expected $(Vs[i]), but got $(dec[i])")
        end
    end
    if decoder.metrics["inactivations"] != K-1
        error("expected $(K-1) inactivations, but got $(decoder.metrics["inactivations"])")
    end
    true
end
for K in [1, 10, 100, 200, 250, 254]
    @test test_dense(K)
end
