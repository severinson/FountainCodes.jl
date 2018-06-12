export LT, LTQ

"LT code parameters."
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

"q-ary LT code parameters."
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
Base.repr{CT,DT}(p::LTQ{CT,DT}) = "LTQ{$CT,DT}($(p.K), $(repr(p.dd)))"

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
function coefficient{CT}(X::Int, j::Int, p::LTQ{CT})
    coef = rand(CT)
    while iszero(coef)
        coef = rand(CT)
    end
    return coef
end

function coefficient{CT<:Float64}(X::Int, i::Int, p::LTQ{CT})
    return randn(CT)/1e10+1
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

doc"generate an LT symbol from the intermediate symbols."
function ltgenerate(C::Vector, X::Int, p::LT)
    d, a, b = trip(X, p)
    while (b >= p.L)
        b = (b + a) % p.Lp
    end
    neighbours = Vector{Int}(min(d, p.L))
    neighbours[1] = b+1
    value = copy(C[b+1])
    for j in 1:min(d-1, p.L-1)
        b = (b + a) % p.Lp
        while (b >= p.L)
            b = (b + a) % p.Lp
        end
        neighbours[j+1] = b+1
        value += C[b+1]
    end
    return BSymbol(X, value, neighbours)
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

"generate an LT symbol from the intermediate symbols."
function ltgenerate{CT}(C::Vector, X::Int, p::LTQ{CT})
    d, _, _ = trip(X, p)
    d = min(d, p.L)
    set = Set{Int}()
    while length(set) < d
        push!(set, rand(1:p.L))
    end
    indices = sort!(collect(set))
    coefficients = [coefficient(X, i, p) for i in 1:d]
    value = sum(C[indices[i]]*coefficients[i] for i in 1:d)
    return QSymbol(X, value, indices, coefficients)
end

"""
    Decoder(p::LT)

Return a decoder for binary LT codes.

"""
function Decoder(p::LT)
    num_buckets = max(3, Int(round(log(p.K))))
    selector = HeapSelect(num_buckets)
    return Decoder{GF256,Vector{GF256},LT,HeapSelect}(
        p,
        selector,
        p.K,
    )
end

"""
    Decoder{CT,DT}(p::LTQ{CT,DT})

Return a decoder for non-binary LT codes.

"""
function Decoder{CT,DT}(p::LTQ{CT,DT})
    num_buckets = max(3, Int(round(log(p.K))))
    selector = HeapSelect(num_buckets)
    return Decoder{CT,Vector{CT},LTQ,HeapSelect}(
        p,
        selector,
        p.K,
    )
end
