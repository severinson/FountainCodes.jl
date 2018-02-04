using Primes

export R10Parameters

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
