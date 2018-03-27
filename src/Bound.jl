# Fountain code bounds

doc"compute log(n!/k!), exactly or approximately, depending on n."
function logfactorial(n::Int, k::Int=1)
    if n < 100
        return logfactorial_exact(n, k)
    else
        return logfactorial_approx(n, k)
    end
end

doc"compute log(n!/k!) exactly."
function logfactorial_exact(n::Int, k::Int=1)
    if n < k
        error("k must be <= n")
    end
    if n < 0 || k < 0
        error("n, k must be non-negative integers")
    end
    r = zero(n)
    for i in k+1:n
        r += log(i)
    end
    return r
end

doc"compute log(n!) approximately. method from [Batir2010]."
function logfactorial_approx(n::Int, k::Int=1)
    if n < k
        error("k must be <= n")
    end
    if n < 0 || k < 0
        error("n, k must be non-negative integers")
    end
    if iszero(n)
        return zero(n)
    end
    r = 1/2 * log(2pi)
    r += n * log(n) - n
    r += 1/2 * log(n+1/6+1/(72n) - 31/(6480n^2) - 139/(155520n^3) + 9871 / (6531840n^4))
    if k != one(k)
        r -= logfactorial_approx(k)
    end
    return r
end

doc"compute log(binomial(n, k))"
function logbinomial(n::Int, k::Int)
    if n < k
        error("k must be <= n")
    end
    if k > (n - k)
        return logfactorial(n, n-k) - logfactorial(k)
    else
        return logfactorial(n, k) - logfactorial(n-k)
    end
end

doc"inner term of ltfailure_lower_reference"
function ltfailure_lower_reference_inner(i::Int, k::Int, epsilon::Number, Omega::Distribution{Univariate, Discrete})
    r = zero(k)
    for d in 1:k
        r += pdf(Omega, d) * binomial(k-i, d) / binomial(k, d)
    end
    r = r^(k*(1+eps))
    return r
end

doc"reference implementation of ltfailure_lower. use only for testing ltfailure_lower."
function ltfailure_lower_reference(k::Int, epsilon::Number, Omega::Distribution{Univariate, Discrete})
    r = zero(k)
    for i in 1:k
        r += (-1)^(i+1)*binom(k, i) * ltfailure_lower_reference_inner(
            i, k, epsilon, Omega,
        )
    end
    return r
end

doc"lower bound on the decoding failure probability of LT codes under ML decoding.
    bound due to [Schotsch2013].
    k - number of input symbols.
    epsilon - relative reception overhead, i.e., 0.2 for a 20% overhead.
    Omega - degree distribution."
function ltfailure_lower(k::Int, epsilon::Number, Omega::Distribution{Univariate, Discrete})

# def nchoosek_log(n, k):
#     '''compute the logarithm of n choose k.'''
#     if k > n:
#         raise ValueError('k must be <= n')
#     if k > (n - k):
#         return factorial_log(n, k=n-k+1) - factorial_log(k)
#     else:
#         return factorial_log(n, k=k+1) - factorial_log(n-k)


