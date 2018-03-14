using Primes, Nulls

export RQParameters

doc"R10 parameters container."
struct RQ <: RaptorCode{NonBinary}
    K::Int # number of source symbols
    Kp::Int # number of source symbols including padding
    S::Int # number of LDPC symbols
    H::Int # number of HDPC symbols
    Hp:: Integer # hamming weight of HDPC symbols
    U::Int # number of PI symbols that are not HDPC symbols
    W::Int # number of LT symbols
    P::Int #number of PI symbols
    L::Integer # =Kp+S+H
    Lp::Integer # used by r10_trip
    function R10Parameters(K::Integer)
        error("not implemented")
    end
end

Base.repr(p::RQ) = "RQ($(p.K))"

# Rand[y, i, m]
# Deg[v]
# Enc[K', C ,(d, a, b, d1, a1, b1)]
# Tuple[Kp, X]

doc"RQ standardized rand function."
function rq_rand(y::Int, i::Int, m::Int) :: Int
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

doc"RQ degree distribution. maps an integer < 2^20 to a degree."
function rq_deg(v::Int) :: Int
    d = Vector(0:30)
    f = [0, 5243, 529531, 704294, 791675, 844104, 879057, 904023, 922747, 937311,
         948962, 958494, 966438, 973160, 978921, 983914, 988283, 992138, 995565,
         998631, 1001391, 1003887, 1006157, 1008229, 1010129, 1011876, 1013490,
         1014983, 1016370, 1017662, 1048576]
    if !(0 <= v <= 1048576)
        error("v must be at least 0 and at most 1048576.")
    end
    j = searchsortedlast(f, v)
    return d[j+1]
end
