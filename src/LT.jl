export LT, LTQ, get_constraint, get_value

"""
    LT{T <: Sampleable{Univariate, Discrete}} <: BinaryCode

Binary Luby Transform erasure correction code. Create a code object
with LT(K, dd).

Arguments

* K::Int: Number of source symbols.

* dd::Sampleable{Univariate, Discrete}: Degree distribution.

"""
struct LT{T <: Sampleable{Univariate, Discrete}} <: BinaryCode
    K::Int # number of source symbols
    L::Int # number of intermediate symbols
    Lp::Int # smallest prime larger than K
    dd::T # degree distribution
    function LT{T}(K::Int, dd::T) where T
        Lp = Primes.nextprime(K)
        new(K, K, Lp, dd)
    end
end
LT(K::Int, dd::T) where {T <: Sampleable{Univariate, Discrete}} = LT{T}(K, dd)
Base.repr(p::LT) = "LT($(p.K), $(repr(p.dd)))"

"""
    LTQ{CT,DT <: Sampleable{Univariate, Discrete}} <: NonBinaryCode

Non-binary Luby Transform erasure correction code with coefficients of
type CT. Create a code object with LT{CT}(K, dd).

Arguments

* K::Int: Number of source symbols.

* dd::Sampleable{Univariate, Discrete}: Degree distribution.

"""
struct LTQ{CT,DT <: Sampleable{Univariate, Discrete}} <: NonBinaryCode
    K::Int # number of source symbols
    L::Int # number of intermediate symbols
    Lp::Int # smallest prime larger than K
    dd::DT # degree distribution
    function LTQ{CT,DT}(K::Int, dd::DT) where {CT,DT}
        Lp = Primes.nextprime(K)
        new(K, K, Lp, dd)
    end
end
function LTQ(K::Int, dd::DT) where DT <: Sampleable{Univariate, Discrete}
    LTQ{GF256,DT}(K, dd)
end
function LTQ{CT}(K::Int, dd::DT) where {CT,DT <: Sampleable{Univariate, Discrete}}
    LTQ{CT,DT}(K, dd)
end
Base.repr(p::LTQ{CT,DT}) where CT where DT = "LTQ{$CT,$DT}($(p.K), $(repr(p.dd)))"

"LT codes have no pre-code, so do nothing."
function precode!(C::Vector, p::Code)
    return C
end

"Map a number 0 <= v <= 1 to a degree."
function deg(v::Real, p::Code) :: Int
    return quantile(p.dd, v)
end

"""
    coefficient{CT,DT}(p::LTQ{CT})

Return a randomly generated coefficient.

"""
function coefficient(X::Int, j::Int, p::LTQ{CT}) where CT
    coef = rand(CT)
    while iszero(coef) coef = rand(CT) end
    return coef
end

"""
    coefficient(X::Int, i::Int, p::LTQ{Float64})

Return a coefficient drawn from a standard Gaussian distribution. See
"Numerically Stable Real Number Codes Based on Random Matrices" by
Chen et al. for a motivation.

"""
function coefficient(X::Int, i::Int, p::LTQ{Float64})
    return randn()
end

"Maps an encoding symbol ID X to a triple (d, a, b)"
function trip(X::Int, p::Union{LT,LTQ})
    Q = 65521 # the largest prime smaller than 2^16
    JK = 1 # no systematic indices for LT codes
    A = (53591 + JK*997) % Q
    B = 10267*(JK+1) % Q
    Y = (B + X*A) % Q
    v = r10_rand(Y, 0, 2<<19) / (2<<19)
    d = deg(v, p)
    a = 1 + r10_rand(Y, 1, p.Lp-1)
    b = r10_rand(Y, 2, p.Lp)
    return d, a, b
end

"""
    get_constraint(lt::LT, X::Integer)

Return the SparseVector corresponding to the X-th LT code constraint.

"""
function get_constraint(lt::LT, X::Integer)
    d, a, b = trip(X, lt)
    while (b >= lt.L)
        b = (b + a) % lt.Lp
    end
    Is = zeros(Int, min(d, lt.L))
    Is[1] = b+1
    for j in 1:min(d-1, lt.L-1)
        b = (b + a) % lt.Lp
        while (b >= lt.L)
            b = (b + a) % lt.Lp
        end
        Is[j+1] = b+1
    end
    Vs = ones(Bool, d)
    return sparsevec(Is, Vs, lt.L)
end

"""
    get_constraint(code::LTQ, X::Integer)

Return the SparseVector corresponding to the X-th LT code constraint.

"""
function get_constraint(ltq::LTQ, X::Integer)
    Random.seed!(X)
    d, _, _ = trip(X, ltq)
    d = min(d, ltq.L)
    set = Set{Int}()
    while length(set) < d
        push!(set, rand(1:ltq.L))
    end
    Is = sort!(collect(set))
    Vs = [coefficient(X, i, ltq) for i in 1:d]
    return sparsevec(Is, Vs, ltq.L)
end

"""
    get_value(code, X, src)

Return the value of the X-th constraint.

"""
function get_value(code, X::Integer, src)
    @assert X >= 0
    return get_value(get_constraint(code, X), src)
end

"""
    get_value(constraint, src)

Return the value of the X-th constraint.

"""
function get_value(constraint::SparseVector, src)
    @assert length(constraint) == length(src)
    @assert length(src) > 0
    rv = zero(src[1])
    for (i, v) in zip(constraint.nzind, constraint.nzval)
        rv += v.*src[i]
    end
    return rv
end

# Using this method results in significantly higher
# probability of decoding failure.
#
# "generate an LT symbol from the intermediate symbols."
# function ltgenerate{CT}(C::Vector, X::Int, p::LTQ{CT})
#     d, a, b = trip(X, p)
#     while (b >= p.L)
#         b = (b + a) % p.Lp
#     end
#     indices = Vector{Int}(min(d, p.L))
#     coefficients = Vector{CT}(min(d, p.L))
#     indices[1] = b+1
#     coefficients[1] = coefficient(X, p)
#     value = C[b+1] * coefficients[1]
#     for j in 1:min(d-1, p.L-1)
#         b = (b + a) % p.Lp
#         while (b >= p.L)
#             b = (b + a) % p.Lp
#         end
#         indices[j+1] = b+1
#         coefficients[j+1] = coefficient(X, p)
#         value = value + C[b+1] * coefficients[j+1]
#     end
#     return QSymbol(X, value, indices, coefficients)
# end

"""
    Decoder{VT}(p::LT)

Return a decoder for binary LT codes with value type VT.

"""
function Decoder(p::LT) where VT
    num_buckets = max(3, Int(round(log(p.K))))
    selector = HeapSelect(num_buckets, p.L)
    return Decoder{Bool}(p, selector, p.K)
end

"""
    Decoder{VT}(p::LTQ{CT}) where {CT,VT}

Return a decoder for non-binary LT codes with coefficient type CT and
value type VT. Note that it must be possible to multiply instances of
VT by instances of CT.

"""
function Decoder(p::LTQ{CT}) where CT
    num_buckets = max(3, Int(round(log(p.K))))
    selector = HeapSelect(num_buckets, p.L)
    return Decoder{CT}(p, selector, p.K)
end
