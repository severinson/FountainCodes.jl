using Primes

export LTParameters

doc"LT parameters."
struct LTParameters <: Parameters
    K::Integer # number of source symbols
    L::Integer # number of intermediate symbols
    Lp::Integer
    dd::Soliton # degree distribution object
    function LTParameters(K::Integer, dd::Soliton)
        if K != dd.K
            error("K = $k != dd.K = $(dd.k)")
        end
        Lp = Primes.nextprime(K)
        new(K, K, Lp, dd)
    end
end

Base.repr(p::LTParameters) = "LTParameters($(p.K), $(repr(p.dd)))"

doc"Map a number 0 <= v <= 1 to a degree."
function deg(p::LTParameters) :: Int
    return rand(p.dd)
end

doc"Maps an encoding symbol ID X to a triple (d, a, b)"
function trip(X::Int, p::LTParameters)
    Q = 65521 # the largest prime smaller than 2^16
    JK = J[p.K+1]
    A = (53591 + JK*997) % Q
    B = 10267*(JK+1) % Q
    Y = (B + X*A) % Q
    d = deg(p)
    a = 1 + r10_rand(Y, 1, p.Lp-1)
    b = r10_rand(Y, 2, p.Lp)
    return d, a, b
end
