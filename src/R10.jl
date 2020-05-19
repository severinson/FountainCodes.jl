export R10, precode

"""
    R10

Raptor10 rateless erasure code. Construct with R10(K::Integer), where
K is the number of source symbols satisfying 4 <= K <= 8192.

"""
struct R10 <: AbstractErasureCode
    K::Int # number of source symbols
    S::Int # number of LDPC symbols
    H::Int # number of HDPC symbols
    Hp:: Int # hamming weight of HDPC symbols
    L::Int # =K+S+H
    Lp::Int # used by r10_trip
    precode::Vector{SparseVector{Bool,Int}} # cached precode constraints
    function R10(K::Integer)
        4 <= K <= 8192 || throw(DomainError(K, "R10 codes support 4 <= K <= 8192."))

        # the smallest positive integer such that X(X-1) >= 2K
        X = ceil(1/2 + sqrt(1/4+2K))

        # S is the smallest prime integer such that S >= ceil(0.01K) + X
        S = Primes.nextprime(Int(ceil(0.01K)+X))

        # H is the smallest integer such that choose(H,ceil(H/2)) >= K + S
        H = 0
        while binomial(H, Int(ceil(H/2))) < (K + S)
            H += 1
        end
        Hp = ceil(H/2)
        L = K+S+H
        Lp = Primes.nextprime(L)
        code = new(K, S, H, Hp, L, Lp, Vector{SparseVector{Bool,Int}}())

        # pre-compute the precode constraints
        ldpc_constraints!(code)
        hdpc_constraints!(code)
        return code
    end
end
Base.show(io::IO, c::R10) = print(io, "R10($(c.K))")
dimension(r10::R10) = r10.L

"""

R10 standardized random number generator. Maps the integers X and i to
a psuedo-random number between 0 and m-1.

"""
function r10_rand(X::Integer, i::Integer, m::Integer)
    X >= 0 || throw(DomainError(X, "ESI must be non-negative."))
    i >= 0 || throw(DomainError(i, "i must be non-negative."))
    m > 0 || throw(DomainError(m, "m must be positive."))
    return UInt32(xor(V0[(X+i) % 256 + 1], V1[(div(X, 256) + i) % 256 + 1]) % m)
end

"""Map a uniformly distributed random number v to a degree."""
function r10_degree(v::Integer)
    0 <= v < 1048576 || throw(DomainError(v, "v must satisfy 0 <= v < 1048576."))
    d = [1, 2, 3, 4, 10, 11, 40]
    f = [0, 10241, 491582, 712794, 831695, 948446, 1032189, 1048576]
    j = 1
    while !(f[j-1+1] <= v < f[j+1]) j += 1 end
    return d[j]
end

"""Map an ESI X to a triple (d, a, b)."""
function trip(X::Integer, code::R10)
    Q = 65521 # the largest prime smaller than 2^16
    J = R10_SYSTEMATIC_INDICES[code.K-4+1] # the first element is for K = 4
    A = (53591 + J*997) % Q
    B = 10267*(J+1) % Q
    Y = (B + X*A) % Q
    v = r10_rand(Y, 0, 2<<19)
    d = r10_degree(v)
    a = 1 + r10_rand(Y, 1, code.Lp-1)
    b = r10_rand(Y, 2, code.Lp)
    return d, a, b
end

"""add the LDPC constraints to the code"""
function ldpc_constraints!(code::R10)

    # the i-th element of indices is the vector of non-zero column
    # indices of the i-th LDPC constraint.
    indices = [Vector{Int}() for _ in 1:code.S]

    # each source symbol participates in exactly 3 LDPC constraints.
    for i in 0:code.K-1
        a = 1 + floor(Int, i/code.S) % (code.S-1)
        b = i % code.S
        push!(indices[b+1], i+1)
        b = (b + a) % code.S
        push!(indices[b+1], i+1)
        b = (b + a) % code.S
        push!(indices[b+1], i+1)
    end

    # at this point the i-th constraint has non-zero value. adding the
    # index of the constraint makes it a parity check constraint,
    # i.e., it's guaranteed to have value zero.
    for (i, Is) in zip((code.K+1):(code.K+code.S), indices)
        push!(Is, i)
        sort!(Is)

        # add to precode constraints
        push!(code.precode, sparsevec(Is, ones(Bool, length(Is)), code.L))
    end
    return
end

"""add the HDPC constraints to the code"""
function hdpc_constraints!(code::R10)
    for h in 0:code.H-1 # loop over HDPC rows
        g = 0
        Is = Vector{Int}() # non-zero indices of h-th HDPC constraint
        for j in 0:code.K+code.S-1 # loop over source/LDPC symbols
            g = nextgray(g, code.Hp)
            if !iszero(g & (1 << h)) push!(Is, j+1) end
        end
        push!(Is, code.K+code.S+h+1)
        sort!(Is)
        push!(code.precode, sparsevec(Is, ones(Bool, length(Is)), code.L))
    end
    return
end

"""
    get_constraint(p::R10, X::Integer)

Return the SparseVector corresponding to the X-th LT code constraint.

"""
function get_constraint(p::R10, X::Integer)

    # precode constraints use negative ESIs
    if X < 0 return p.precode[-1*X] end

    # compute LT constraint with ESI X
    d, a, b = trip(X, p)
    while (b >= p.L)
        b = (b + a) % p.Lp
    end
    Is = zeros(Int, min(d, p.L))
    Is[1] = b+1
    for j in 1:min(d-1, p.L-1)
        b = (b + a) % p.Lp
        while (b >= p.L)
            b = (b + a) % p.Lp
        end
        Is[j+1] = b+1
    end
    Vs = ones(Bool, min(d, p.L))
    return sparsevec(Is, Vs, p.L)
end

"""
    precode!(src, code::R10)

Compute the intermediate symbols and store them in src.

"""
function precode(src::Vector{VT}, code::R10) where VT
    length(src) == code.K || throw(DimensionMismatch("src must have length $(code.K)"))
    Vs = [zero(src[1]) for _ in 1:(code.S+code.H)] # parity symbol values
    append!(Vs, src) # source symbol values
    return decode(code, -(code.S+code.H):code.K-1, Vs)
end

"""
    Decoder(code::R10)

Return a decoder for Raptor10 codes.

"""
function Decoder(code::R10)
    num_buckets = 41 # 1 more than the max LT degree
    selector = HeapSelect(num_buckets, code.L)
    d = Decoder{Bool}(selector, code.L)
    return d
end
