using Primes, Distributions

export LTParameters, QLTParameters

doc"LT code parameters."
struct LTParameters{T <: Sampleable{Univariate, Discrete}} <: LTCode{Binary}
    K::Integer # number of source symbols
    L::Integer # number of intermediate symbols
    Lp::Integer
    dd::T # degree distribution
    function LTParameters{T}(K::Integer, dd::T) where T
        Lp = Primes.nextprime(K)
        new(K, K, Lp, dd)
    end
end

LTParameters(K::Integer, dd::T) where {T <: Sampleable{Univariate, Discrete}} = LTParameters{T}(K, dd)

Base.repr(p::LTParameters) = "LTParameters($(p.K), $(repr(p.dd)))"

doc"q-ary LT code parameters."
struct QLTParameters{DT <: Sampleable{Univariate, Discrete}, CT} <: LTCode{NonBinary}
    K::Integer # number of source symbols
    L::Integer # number of intermediate symbols
    Lp::Integer
    dd::DT # degree distribution
    function QLTParameters{DT,CT}(K::Integer, dd::DT) where {DT, CT}
        Lp = Primes.nextprime(K)
        new(K, K, Lp, dd)
    end
end

function QLTParameters(K::Integer, dd::DT) where DT <: Sampleable{Univariate, Discrete}
    QLTParameters{DT,GF256}(K, dd)
end

doc"Map a number 0 <= v <= 1 to a degree."
function deg(v::Real, p::LTCode) :: Int
    return quantile(p.dd, v)
end

doc"Map a number 0 <= v <= 1 to a coefficient."
function coefficient{DT,CT}(p::QLTParameters{DT,CT})
    return rand(CT)
    # return one(CT)
end

doc"Maps an encoding symbol ID X to a triple (d, a, b)"
function trip(X::Int, p::LTCode)
    Q = 65521 # the largest prime smaller than 2^16
    JK = J[p.K+1]
    A = (53591 + JK*997) % Q
    B = 10267*(JK+1) % Q
    Y = (B + X*A) % Q
    v = r10_rand(Y, 0, 2<<19) / (2<<19)
    d = deg(v, p)
    a = 1 + r10_rand(Y, 1, p.Lp-1)
    b = r10_rand(Y, 2, p.Lp)
    return d, a, b
end

doc"Generate an LT symbol from the intermediate symbols."
function lt_generate(C::Vector, X::Int, p::LTCode{Binary})
    d, a, b = trip(X, p)
    while (b >= p.L)
        b = (b + a) % p.Lp
    end
    neighbours = Vector{Int}(min(d, p.L))
    neighbours[1] = b+1
    value = C[b+1]
    for j in 1:min(d-1, p.L-1)
        b = (b + a) % p.Lp
        while (b >= p.L)
            b = (b + a) % p.Lp
        end
        neighbours[j+1] = b+1
        value = value + C[b+1]
    end
    return R10Symbol(X, value, neighbours)
end

doc"Generate an LT symbol from the intermediate symbols."
function lt_generate(C::Vector, X::Int, p::LTCode{NonBinary})
    d, a, b = trip(X, p)
    while (b >= p.L)
        b = (b + a) % p.Lp
    end
    indices = Vector{Int}(min(d, p.L))
    coefficients = Vector{GF256}(min(d, p.L))
    indices[1] = b+1
    coefficients[1] = coefficient(p)
    value = C[b+1] * coefficients[1]
    for j in 1:min(d-1, p.L-1)
        b = (b + a) % p.Lp
        while (b >= p.L)
            b = (b + a) % p.Lp
        end
        indices[j+1] = b+1
        coefficients[j+1] = coefficient(p)
        value = value + C[b+1] * coefficients[j+1]
    end
    return QSymbol(X, value, indices, coefficients)
end
