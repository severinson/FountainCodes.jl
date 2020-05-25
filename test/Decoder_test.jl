# Basic decoder tests

using FountainCodes, LinearAlgebra, SparseArrays, Test

"""Test decoding with a diagonal constraint matrix."""
function test_diagonal(CT, K::Integer)
    constraints = [sparsevec([i], [one(CT)], K) for i in 1:K]
    Vs = [GF256(i%256) for i in 0:K-1]
    d = Decoder{CT}(K)
    dec = decode(constraints, Vs, decoder=d)
    for i in 1:K
        if dec[i] != Vs[i]
            error("expected $(Vs[i]), but got $(dec[i])")
        end
    end
    if d.metrics["inactivations"] != 0
        error("expected 0 inactivations, but got $(d.metrics["inactivations"])")
    end
    return true
end
for K in [1, 10, 100, 200, 250, 254]
    @test test_diagonal(Bool, K)
end
for K in [1, 10, 100, 200, 250, 254]
    @test test_diagonal(GF256, K)
end

"""Test decoding for a dense constraint matrix."""
function test_dense(K::Integer)
    Is = collect(1:K)
    Cs = [GF256(i) for i in 1:K]
    constraints = [sparsevec(Is, circshift(Cs, i), K) for i in 0:K-1]
    src = [GF256(i%256) for i in 0:K-1]
    enc = [dot(constraint, src) for constraint in constraints]
    d = Decoder{GF256}(K)
    dec = decode(constraints, enc, decoder=d)
    for i in 1:K
        if dec[i] != src[i]
            error("expected $(Vs[i]), but got $(dec[i])")
        end
    end
    if d.metrics["inactivations"] != K-1
        error("expected $(K-1) inactivations, but got $(d.metrics["inactivations"])")
    end
    return true
end
for K in [1, 10, 100, 200, 250, 254]
    @test test_dense(K)
end
