export R10Parameters, precode!, ltgenerate

doc"R10 parameters container."
struct R10Parameters <: RaptorCode{Binary}
    K::Integer # number of source symbols
    S::Integer # number of LDPC symbols
    H::Integer # number of HDPC symbols
    Hp:: Integer # hamming weight of HDPC symbols
    L::Integer # =K+S+H
    Lp::Integer # used by r10_trip
    function R10Parameters(K::Integer)
        if !(4 <= K <= 8192)
            error("R10 codes support 4 <= K <= 8192 source symbols.")
        end

        # X is the smallest positive integer such that X(X-1) >= 2K
        X = ceil(1/2 + sqrt(1/4+2K))

        # S is the smallest prime integer such that S >= ceil(0.01K) + X
        S = Primes.nextprime(Int64(ceil(0.01K)+X))

        # H is the smallest integer such that choose(H,ceil(H/2)) >= K + S
        H = 0
        while binomial(H,Int64(ceil(H/2))) < K + S
            H += 1
        end
        Hp = ceil(H/2)
        L = K+S+H
        Lp = Primes.nextprime(L)
        new(K, S, H, Hp, L, Lp)
    end
end

Base.repr(p::R10Parameters) = "R10Parameters($(p.K))"

doc"pre-code input data. C must be an array of source symbols of length L."
function precode!(C::Vector, p::R10Parameters)
    C = r10_ldpc_encode!(C, p, null)
    C = r10_hdpc_encode!(C, p, null)
    return C
end

doc"pre-code input data. C must be an array of source symbols of length L."
function precode!(C::Vector, p::R10Parameters, N::Vector)
    C = r10_ldpc_encode!(C, p, N)
    C = r10_hdpc_encode!(C, p, N)
    return C
end

doc"generate an LT symbol from the intermediate symbols."
function ltgenerate(C::Vector, X::Int, p::Code)
    d, a, b = trip(X, p)
    while (b >= p.L)
        b = (b + a) % p.Lp
    end
    N = Vector{Int}(min(d, p.L))
    N[1] = b+1
    value = C[b+1]
    for j in 1:min(d-1, p.L-1)
        b = (b + a) % p.Lp
        while (b >= p.L)
            b = (b + a) % p.Lp
        end
        N[j+1] = b+1
        value = value + C[b+1]
    end
    return BSymbol(X, value, N)
end

doc"Generate R10 precode LDPC symbols in-place at N (K+1) to (K+S)."
function r10_ldpc_encode!(C::Vector, p::R10Parameters, N=null)
    if length(C) != p.L
        error("C must have length p.L = $p.L")
    end
    if !(N isa Null) && length(N) != p.L
        error("N must have length p.L = $p.L")
    end
    for i in 1:p.S
        C[p.K+i] = zero(C[1])
    end
    for i in 0:p.K-1
        v = C[i+1]
        a = 1 + Int64((floor(i/p.S) % (p.S-1)))
        b = i % p.S
        C[p.K+b+1] = C[p.K+b+1] + v
        if !(N isa Null)
            push!(N[p.K+b+1], i+1)
        end
        b = (b + a) % p.S
        C[p.K+b+1] = C[p.K+b+1] + v
        if !(N isa Null)
            push!(N[p.K+b+1], i+1)
        end
        b = (b + a) % p.S
        C[p.K+b+1] = C[p.K+b+1] + v
        if !(N isa Null)
            push!(N[p.K+b+1], i+1)
        end
    end
    return C
end

doc"Generate R10 precode HDPC symbols in-place at N (K+S+1) to (K+S+H)."
function r10_hdpc_encode!(C::Vector, p::R10Parameters, N=null)
    if length(C) != p.L
        error("C must have length p.L = $p.L")
    end
    if !(N isa Null) && length(N) != p.L
        error("N must have length p.L = $p.L")
    end
    for i in 1:p.H
        C[p.K+p.S+i] = zero(C[1])
    end
    for h in 0:p.H-1
        j = 0
        for g in gray(p.K+p.S, p.Hp)
            if !iszero(g & (1 << h))
                C[p.K+p.S+h+1] = C[p.K+p.S+h+1] + C[j+1]
                if !(N isa Null)
                    push!(N[p.K+p.S+h+1], j+1)
                end
            end
            j += 1
        end
    end
    return C
end

doc"R10 standardized rand function."
function r10_rand(X::Int, i::Int, m::Int) :: Int
    if X < 0
        error("X must be non-negative")
    end
    if i < 0
        error("i must be non-negative")
    end
    if m <= 0
        error("m must be positive")
    end
    return xor(
        V0[(X+i) % 256 + 1],
        V1[(div(X, 256) + i) % 256 + 1]
    ) % m
end

doc"Map a uniformly distributed random number v to a degree."
function deg(v::Int) :: Int
    d = [1, 2, 3, 4, 10, 11, 40]
    f = [10241, 491582, 712794, 831695, 948446, 1032189, 1048576]
    if !(0 <= v <= 1048576)
        error("v must be at least 0 and at most 1048576.")
    end
    j = 1
    while v > f[j]
        j += 1
    end
    return d[j]
end

doc"Maps an encoding symbol ID X to a triple (d, a, b)"
function trip(X::Int, p::R10Parameters)
    Q = 65521 # the largest prime smaller than 2^16
    JK = J[p.K+1]
    A = (53591 + JK*997) % Q
    B = 10267*(JK+1) % Q
    Y = (B + X*A) % Q
    v = r10_rand(Y, 0, 2<<19)
    d = deg(v)
    a = 1 + r10_rand(Y, 1, p.Lp-1)
    b = r10_rand(Y, 2, p.Lp)
    return d, a, b
end
