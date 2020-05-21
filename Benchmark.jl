using FountainCodes

function benchmark_r10(K=1000, r=1200, m=256, n=100)
    r10 = R10(K)
    d = Decoder(r10)
    src = [GF256(i % 256) for i in 1:K]
    inter = precode(src, r10)
    Vs = [zero(inter[1]) for _ in 1:(r10.S+r10.H)] # parity symbol values
    append!(Vs, [get_value(r10, X, inter) for X in 0:r-1]) # LT symbol values
    Xs = -(r10.S+r10.H):r-1 # ESIs incl. parity symbols

    t = 0.0
    for _ in 1:n
        t += @elapsed decode(r10, Xs, Vs)
    end
    t /= n
    println(r10)
    println("Decoding time: $t s")
end

function benchmark_lt(K=1000, r=1300, m=256, n=100, M=K-1, δ=1e-6)
    dd = Soliton(K, M, δ)
    lt = LT(K, dd)
    src = [Vector{GF256}([i % 256]) for i in 1:K]
    Vs = [get_value(lt, X, src) for X in 1:r]
    t = 0.0
    for _ in 1:n
        t += @elapsed decode(lt, 1:r, Vs)
    end
    t /= n
    println(lt)
    println("Decoding time: $t s")
end