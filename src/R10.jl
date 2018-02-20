using Primes

export R10Parameters

# R10 encoding instructions:
# - Allocate an array C of length L.
# - Assign the K source symbols to the first K entries of C.
# - Call r10_ldpc_encode with C as its argument to generate the LDPC symbols.
# - Call r10_hdpc_encode with C as its argument to generate the HDPC symbols.
# - Call r10_lt_encode with C as its argument to generate an LT symbol.

doc"R10 parameters container."
struct R10Parameters <: Parameters
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

doc"Generate R10 precode LDPC symbols in-place at indices (K+1) to (K+S)."
function r10_ldpc_encode!{VT<:Value}(C::Vector{ISymbol{VT}}, p::R10Parameters)
    if length(C) != p.L
        error("C must have length p.L = $p.L")
    end
    neighbours = [Set{Int}() for _ in 1:p.S]
    values = [R10Value(0) for _ in 1:p.S]
    for i in 0:p.K-1
        v = C[i+1].value
        a = 1 + Int64((floor(i/p.S) % (p.S-1)))
        b = i % p.S
        push!(neighbours[b+1], i+1)
        values[b+1] = values[b+1] + v
        b = (b + a) % p.S
        push!(neighbours[b+1], i+1)
        values[b+1] = values[b+1] + v
        b = (b + a) % p.S
        push!(neighbours[b+1], i+1)
        values[b+1] = values[b+1] + v
    end
    for i in 1:p.S
        C[p.K+i] = ISymbol(values[i], neighbours[i])
    end
end

doc"Generate R10 precode HDPC symbols in-place at indices (K+S+1) to (K+S+H)."
function r10_hdpc_encode!{VT<:Value}(C::Vector{ISymbol{VT}}, p::R10Parameters)
    if length(C) != p.L
        error("C must have length p.L = $p.L")
    end
    neighbours = [Set{Int}() for _ in 1:p.H]
    values = [R10Value(0) for _ in 1:p.H]
    for h in 0:p.H-1
        j = 0
        for g in gray(p.K+p.S, p.Hp)
            if g & (1 << h) != 0
                push!(neighbours[h+1], j+1)
                values[h+1] = values[h+1] + C[j+1].value
            end
            j += 1
        end
    end
    for i in 1:p.H
        C[p.K+p.S+i] = ISymbol(values[i], neighbours[i])
    end
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

doc"Generate an LT symbol from the intermediate symbols."
function lt_generate{VT<:Value}(C::Vector{ISymbol{VT}}, X::Int, p::Parameters)
    d, a, b = trip(X, p)
    while (b >= p.L)
        b = (b + a) % p.Lp
    end
    neighbours = Vector{Int}(min(d, p.L))
    neighbours[1] = b+1
    value = C[b+1].value
    for j in 1:min(d-1, p.L-1)
        b = (b + a) % p.Lp
        while (b >= p.L)
            b = (b + a) % p.Lp
        end
        neighbours[j+1] = b+1
        value = value + C[b+1].value
    end
    return R10Symbol(X, value, neighbours)
end
