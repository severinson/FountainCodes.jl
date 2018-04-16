export RQ, precode!, ltgenerate

"""
    RQ

RaptorQ code.

"""
struct RQ <: NonBinaryCode
    K::Int # number of source symbols
    Kp::Int # number of source symbols including padding
    J::Int # systematic index
    S::Int # number of LDPC symbols
    H::Int # number of HDPC symbols
    W::Int # number of LT symbols
    L::Int # =Kp+S+H
    P::Int # number of PI symbols
    P1::Int # smallest prime >= P
    U::Int # number of PI symbols that are not HDPC symbols
    B::Int # number of LT symbols excluding the LDPC symbols
    function RQ(K::Integer)
        Kp, J, S, H, W = RQ_parameters(K)
        L = Kp + S + H
        P = L - W
        P1 = nextprime(P)
        U = P - H
        B = W - S
        return new(K, Kp, J, S, H, W, L, P, P1, U, B)
    end
end

Base.repr(p::RQ) = "RQ($(p.K))"

"""
    RQ_parameters(K)

Return the parameter tuple (Kp, J, S, H, W) for RaptorQ codes with K source
symbols.

"""
function RQ_parameters(K::Int)
    if K > 56403
        error("there can be at most 56403 source symbols")
    end
    i = start(searchsorted(RQ_parameter_table[:, 1], K))
    Kp = RQ_parameter_table[i, 1]
    J = RQ_parameter_table[i, 2]
    S = RQ_parameter_table[i, 3]
    H = RQ_parameter_table[i, 4]
    W = RQ_parameter_table[i, 5]
    return Kp, J, S, H, W
end

"""
    RQ_rand(y::Int, i::Int, m::Int)

RaptorQ standardized rand function. Maps the integers y and i to a
psuedo-randomly generated number between 0 and m-1.

"""
function RQ_rand(y::Int, i::Int, m::Int) :: Int
    if y < 0
        error("y must be non-negative")
    end
    if !(0 <= i < 256)
        error("i must be non-negative and less than 256")
    end
    if m <= 0
        error("m must be positive")
    end
    x0 = (y + i) % 256
    x1 = (Int(floor(y / 256)) + i) % 256
    x2 = (Int(floor(y / 65536)) + i) % 256
    x3 = (Int(floor(y / 16777216)) + i) % 256
    result = xor(V0[x0+1], V1[x1+1])
    result = xor(result, V2[x2+1])
    result = xor(result, V3[x3+1])
    return result % m
end

"""
    RQ_deg(v::Int)

RaptorQ degree distribution. Maps an integer 0 <= v < 2^20 to a degree. v must
be uniformly distributed over the range 0 to 2^20.

"""
function RQ_deg(v::Int) :: Int
    d = Vector(0:30)
    f = [0, 5243, 529531, 704294, 791675, 844104, 879057, 904023, 922747, 937311,
         948962, 958494, 966438, 973160, 978921, 983914, 988283, 992138, 995565,
         998631, 1001391, 1003887, 1006157, 1008229, 1010129, 1011876, 1013490,
         1014983, 1016370, 1017662, 1048576]
    if !(0 <= v < 1048576)
        error("v must be non-negative and less than 1048576.")
    end
    j = searchsortedlast(f, v)
    return d[j+1]
end

"""
    RQ_tuple(Kp::Int, X::Int)

RaptorQ standardized tuple function. Takes as input the number of source symbols
Kp and a symbol identifier X and returns the tuple (d, a, b, d1, a1, b1). These
numbers uniquely determine the LT symbol with ISI X.

"""
function RQ_tuple(X::Int, c::RQ)
    A = 53591 + 997c.J
    if (A % 2 == 0) A += 1 end
    B = 10267*(c.J+1)
    y = (B + X*A) % 2^32
    v = RQ_rand(y, 0, 2^20)
    d = min(RQ_deg(v), c.W-2)
    a = 1 + RQ_rand(y, 1, c.W-1)
    b = RQ_rand(y, 2, c.W)
    if (d < 4) d1 = 2 + RQ_rand(X, 3, 2)
    else d1 = 2 end
    a1 = 1 + RQ_rand(X, 4, c.P1-1)
    b1 = RQ_rand(X, 5, c.P1)
    return d, a, b, d1, a1, b1
end

"""
    ltgenerate(C::Vector, X::Int, c::RQ)

Generate an LT symbol from the intermediate symbol vector C and its encoded
symbol identifier, X.

"""
function ltgenerate(C::Vector, X::Int, c::RQ)
    d, a, b, d1, a1, b1 = RQ_tuple(X, c)
    indices = Vector{Int}(d+d1)
    value = C[b+1]
    indices[1] = b+1
    for j in 1:d-1
        b = (b + a) % c.W
        value += C[b+1]
        indices[j+1] = b+1
    end
    while (b1 >= c.P) b1 = (b1+a1) % c.P1 end
    value += C[c.W+b1+1]
    indices[d+1] = c.W+b1+1
    for j in 1:d1-1
        b1 = (b1 + a1) % c.P1
        while (b1 >= c.P) b1 = (b1+a1) % c.P1 end
        value += C[c.W+b1+1]
        indices[d+j+1] = c.W+b1+1
    end
    sort!(indices)
    return BSymbol(X, value, indices)
end

"""
    precode!(C::Vector, c::RQ, N=missing)

RaptorQ precode. C is assumed to be a vector of length L, where the elements at
index 1, ..., K are the source symbols. After calling this method elements K+1,
..., K+S will be the LDPC symbols, and elements K+S+1, ..., K+S+H will be the
HDPC symbols.

# TODO: separate relations and value generation for LT symbols since we often
need only the indices and not the value.

"""
function precode!(C::Vector, c::RQ)
    d = Decoder(c)
    if length(C) != c.L
        error("C must have length L")
    end
    for i in c.K+1:c.L
        C[i] = zero(C[1])
    end
    for X in 1:c.Kp
        s = ltgenerate(C, X, c)
        s = BSymbol(X, C[X], s.neighbours)
        add!(d, s)
    end
    decode!(d)
    get_source!(C, d)
    return C
end

"""
    precode_relations(C::Vector, c::RQ, N)

Return a Vector of tuples [(indices, coefficients), ...], describing the
intermediate symbol constraints. The first S entries correspond to LDPC
relations and the remaining H entries correspond to HDPC relations.

"""
function precode_relations(c::RQ)
    N = Vector{Tuple{Vector{Int},Vector{GF256}}}(c.L)
    N = RQ_ldpc_constraints!(N, c)
    N = RQ_hdpc_constraints!(N, c)
    return N
end

function RQ_ldpc_constraints!(N::Vector, c::RQ)
    if length(N) != c.L
        error("N must have length c.L = $(c.L)")
    end
    for i in 1:c.S
        N[i] = (Vector{Int}(), Vector{GF256}())
    end
    for i in 0:c.B-1
        a = 1 + Int(floor(i/c.S))
        b = i % c.S
        push!(N[b+1][1], i+1)
        push!(N[b+1][2], true)
        b = (b + a) % c.S
        push!(N[b+1][1], i+1)
        push!(N[b+1][2], true)
        b = (b + a) % c.S
        push!(N[b+1][1], i+1)
        push!(N[b+1][2], true)
    end
    for i in 0:c.S-1
        a = i % c.P
        b = (i + 1) % c.P
        push!(N[i+1][1], i+c.Kp+1)
        push!(N[i+1][2], true)
        push!(N[i+1][1], c.W+a+1)
        push!(N[i+1][2], true)
        push!(N[i+1][1], c.W+b+1)
        push!(N[i+1][2], true)
    end
    return N
end

function RQ_hdpc_constraints!(N::Vector, c::RQ)
    alpha = 0x02
    if length(N) != c.L
        error("N must have length c.L = $(c.L)")
    end
    i0 = c.S+1
    for i in i0:c.S+c.H
        N[i] = (collect(1:c.Kp+c.S+1), zeros(GF256, c.Kp+c.S+1))
    end
    MT = zeros(GF256, c.H, c.Kp+c.S)
    GAMMA = zeros(GF256, c.Kp+c.S, c.Kp+c.S)
    for j in 0:c.Kp+c.S-2
        i = RQ_rand(j+1, 6, c.H)
        MT[i+1,j+1] = one(GF256)
        i = (i + RQ_rand(j+1, 7, c.H-1) + 1) % c.H
        MT[i+1,j+1] = one(GF256)
    end
    j = c.Kp+c.S
    for i in 1:c.H
        MT[i,j] = alpha^(i-1)
    end
    for i in 0:c.Kp+c.S-1
        for j in 0:i
            GAMMA[i+1,j+1] = alpha^(i-j)
        end
    end
    A::Matrix{GF256} = MT*GAMMA
    for i in 1:c.H
        for j in 1:c.Kp+c.S
            N[i0+i-1][2][j] = A[i,j]
        end
        N[i0+i-1][1][end] = c.S+c.H+i
        N[i0+i-1][2][end] = one(GF256)
    end

    # TODO: more efficient implementation. exploits the structure of the matrix.
    # A = zeros(GF256, c.H, c.Kp+c.S)
    # for j in 0:c.Kp+c.S-2 # col index of MT
    #     i = RQ_rand(j+1, 6, c.H) # row index of MT
    #     l = j # row index of GAMMA
    #     for k in 0:l # col index of GAMMA
    #         A[i+1,k+1] += alpha^(l-k)
    #         N[i0+i][2][k+1] += alpha^(l-k)
    #     end
    #     i = (i + RQ_rand(j+1, 7, c.H-1) + 1) % c.H
    #     for k in 0:l
    #         A[i+1,k+1] += alpha^(l-k)
    #         N[i0+i][2][k+1] += alpha^(l-k)
    #     end
    # end
    # j = c.S+c.H-1
    # l = j
    # for i in 0:c.H-1
    #     for k in 0:i-1
    #         A[i+1,k+1] += alpha^(i+l-k)
    #         N[i0+i][2][k+1] += alpha^(i+l-k)
    #     end
    # end

    return N
end

"""
    Decoder(c::RQ)

Create a RaptorQ decoder and add the relevant constraint symbols.

"""
function Decoder(c::RQ)
    selector = SelectBucket(31)
    d = Decoder{Union{BRow,QRow{GF256}},Vector{GF256},RQ,SelectBucket}(c, selector)
    C = zeros(GF256, c.L)
    N = precode_relations(c)

    # LDPC constraints
    for (indices, _) in view(N, 1:c.S)
        s = BSymbol(-1, Vector{GF256}(), indices)
        add!(d, s)
    end

    # HDPC constrains
    for (indices, coefs) in view(N, c.S+1:c.S+c.H)
        s = QSymbol(-1, Vector{GF256}(), indices, coefs)
        add!(d, s)
    end

    # permanent inactivations
    for i in c.L-c.P+1:c.L
        inactivate!(d, i)
    end

    # only count dynamic inactivations
    push!(d.metrics, "inactivations", -c.P)
    return d
end
