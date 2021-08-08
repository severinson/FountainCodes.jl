using FountainCodes, Test, Distributions

"""return an LT code object and a vector of source symbols"""
function init(K; Tv=GF256, M=K, δ=1e-9)
    dd = Soliton(K, M, δ)
    LT{Tv}(K, dd)
end

function test_decoding(K::Integer, N::Integer; Tv=GF256)
    rng = MersenneTwister(123)
    code = init(K; Tv)
    src = rand_nonzero(rng, Tv, K)
    G = generator_matrix(code, 1:N)
    enc = G'*src
    dec = decode(G, enc)
    @test dec == src
end
test_decoding(100, 120)