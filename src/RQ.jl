export RQ, precode!, ltgenerate

"""

Return the parameter tuple `(Kp, J, S, H, W)` for RaptorQ codes with `K` source symbols.
"""
function RQ_parameters(K::Integer)
    0 < K <= 56403 || throw(DomainError(K, "K must be in [1, 56403]"))
    i = searchsortedfirst(view(RQ_parameter_table, :, 1), K)
    Kp = RQ_parameter_table[i, 1]
    J = RQ_parameter_table[i, 2]
    S = RQ_parameter_table[i, 3]
    H = RQ_parameter_table[i, 4]
    W = RQ_parameter_table[i, 5]
    Kp, J, S, H, W
end

"""
    RQ

RaptorQ code.

"""
struct RQ <: AbstractErasureCode
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
        new(K, Kp, J, S, H, W, L, P, P1, U, B)
    end
end
Base.show(io::IO, code::RQ) = print(io, "RQ($(code.K))")

"""

RaptorQ standardized rand function. Maps the integers `y` and `i` to a psuedo-randomly generated
number between `0` and `m-1` (inclusive).
"""
function RQ_rand(y::Integer, i::Integer, m::Integer)
    0 <= y || throw(DomainError(y, "y must be non-negative"))
    0 <= i < 256 || throw(DomainError(i, "i must be in [0, 256)"))
    0 <= m || throw(DomainError(m, "m must be positive"))
    x0 = (y + i) % 256
    x1 = (Int(floor(y / 256)) + i) % 256
    x2 = (Int(floor(y / 65536)) + i) % 256
    x3 = (Int(floor(y / 16777216)) + i) % 256
    result = xor(V0[x0+1], V1[x1+1])
    result = xor(result, V2[x2+1])
    result = xor(result, V3[x3+1])
    result % m
end

"""

RaptorQ degree distribution. Maps an integer `v` in `[0, 2^20)` to a degree, where `v` is assumed 
to be uniformly distributed over the range.
"""
function RQ_deg(v::Integer)
    0 <= v < 1048576 || throw(DomainError(v, "v must be in [0, 2^20)"))
    searchsortedlast(RQ_DEGREE_TABLE, v)
end

"""

RaptorQ standardized tuple function. Takes as input the number of source symbols `Kp` and an ESI
`X`, and returns the tuple `(d, a, b, d1, a1, b1)`, which uniquely determines the LT symbol with 
ESI `X`.
"""
function RQ_tuple(X::Integer, c::RQ)
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
    d, a, b, d1, a1, b1
end

function ldpc_constraint_matrix(c::RQ)
    indices = [Vector{Int}() for _ in 1:c.S]
    for i in 0:c.B-1
        a = 1 + Int(floor(i/c.S))
        b = i % c.S
        push!(indices[b+1], i+1)
        b = (b + a) % c.S
        push!(indices[b+1], i+1)        
        b = (b + a) % c.S
        push!(indices[b+1], i+1)        
    end
    for i in 0:c.S-1
        a = i % c.P
        b = (i + 1) % c.P
        if c.W+a == i+c.Kp
            push!(indices[i+1], c.W+b+1)
        elseif c.W+b == i+c.Kp
            push!(indices[i+1], c.W+a+1)
        else
            push!(indices[i+1], i+c.Kp+1)
            push!(indices[i+1], c.W+a+1)
            push!(indices[i+1], c.W+b+1)            
        end
    end
    hcat([sparsevec(Is, ones(GF256, length(Is)), c.L) for Is in indices]...)
end

function hdpc_constraint(c::RQ, A::Matrix{GF256}, i::Integer)
    Is = Vector{Int}()
    Vs = Vector{GF256}()
    for j in 1:c.Kp+c.S
        if !iszero(A[i, j])
            push!(Is, j)
            push!(Vs, A[i, j])
        end
    end
    push!(Is, c.Kp+c.S+i)
    push!(Vs, one(GF256))
    sparsevec(Is, Vs, c.L)
end

function hdpc_constraint_matrix(c::RQ)
    alpha = GF256(0x02)
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
    return hcat([hdpc_constraint(c, A, i) for i in 1:c.H]...)

    # TODO: more efficient implementation. exploits the structure of the matrix.
    # (currently doesn't work)
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
end

"""

Return the constraint matrix composed of the LDPC and HDPC constraint.
"""
precode_constraint_matrix(code::RQ) = hcat(ldpc_constraint_matrix(code), hdpc_constraint_matrix(code))

"""

Return a `SparseVector` corresponding to the LT symbol with ESI `X`.
"""
function lt_constraint(c::RQ, X::Integer)
    d, a, b, d1, a1, b1 = RQ_tuple(X, c)
    Is = zeros(Int, d+d1)
    Is[1] = b+1
    for j in 1:d-1
        b = (b + a) % c.W
        Is[j+1] = b+1
    end    
    while b1 >= c.P
        b1 = (b1+a1) % c.P1 
    end    
    Is[d+1] = c.W+b1+1
    for j in 1:d1-1
        b1 = (b1 + a1) % c.P1
        while b1 >= c.P
            b1 = (b1+a1) % c.P1 
        end
        Is[d+j+1] = c.W+b1+1
    end
    sparsevec(Is, ones(GF256, d+d1), c.L)
end

generator_matrix(code::RQ, Xs::AbstractVector{<:Integer}) = hcat([lt_constraint(code, X) for X in Xs]...)
constraint_matrix(code::RQ, Xs::AbstractVector{<:Integer}) = hcat(precode_constraint_matrix(code), generator_matrix(code, Xs))