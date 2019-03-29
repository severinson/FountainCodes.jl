# Fountain code bounds

export ltfailure_lower

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

    # compute last term one step at a time. watch for overflow etc.
    # r += 1/2 * log(n+1/6+1/(72n) - 31/(6480n^2) - 139/(155520n^3) + 9871 / (6531840n^4))    
    r2 = n + 1/6
    tmp = 1/(72n)
    if !isnan(tmp) && !isinf(tmp) && tmp > 0
        r2 += tmp
    end
    tmp = 31/(6480n^2)
    if !isnan(tmp) && !isinf(tmp) && tmp > 0
        r2 -= tmp
    end
    tmp = 139/(155520n^3)
    if !isnan(tmp) && !isinf(tmp) && tmp > 0
        r2 -= tmp
    end    
    tmp = 9871 / (6531840n^4)
    if !isnan(tmp) && !isinf(tmp) && tmp > 0
        r2 += tmp
    end
    r += 1/2 * log(r2)
    
    if k != one(k)
        r -= logfactorial_approx(k)
    end
    return r
end

doc"compute log(binomial(n, k))"
function logbinomial(n::Int, k::Int) :: Float64
    if n < k
        return -Inf
    end
    if k == n
        return 0.0
    end
    r = 0.0
    if k > (n - k)
        r = logfactorial(n, n-k) - logfactorial(k)
    else
        r = logfactorial(n, k) - logfactorial(n-k)
    end
    @assert !isnan(r) && !isinf(r) "logbinomial($n,$k)=$r"
    return r
end

doc"inner term of ltfailure_lower_reference."
function ltfailure_lower_reference_inner(i::Int, k::Int, epsilon::Number, Omega::Distribution{Univariate, Discrete})
    r = zero(k)
    for d in 1:k
        r += pdf(Omega, d) * binomial(k-i, d) / binomial(k, d)
    end
    r = r^(k*(1+epsilon))
    return r
end

doc"reference implementation of ltfailure_lower. use only for testing ltfailure_lower."
function ltfailure_lower_reference(k::Int, epsilon::Number, Omega::Distribution{Univariate, Discrete})
    r = zero(k)
    for i in 1:k
        r += (-1)^(i+1)*binomial(k, i) * ltfailure_lower_reference_inner(
            i, k, epsilon, Omega,
        )
    end
    return r
end

doc"inner term of ltfailure_lower."
function ltfailure_lower_inner(i::Int, k::Int, epsilon::Number, Omega::Distribution{Univariate, Discrete}) :: Float64
    r = zero(Float64)
    for d in 1:k
        if k-i < d || k < d
            continue
        end
        r += exp(log(pdf(Omega, d)) + logbinomial(k-i, d) - logbinomial(k, d))
        if isnan(r)
            println("r=$r")
            println(pdf(Omega,d))
            println(log(pdf(Omega, d)))
            println(logbinomial(k-i, d))
            println(logbinomial(k, d))
            error("foo")
        end
    end
    return log(r) * (k*(1+epsilon))
end

"""
    ltfailure_lower(k::Int, epsilon::Number, Omega::Distribution{Univariate, Discrete})

Lower-bound the decoding failure probability of LT codes with `k` input symbols,
degree distribution `Omega`, at a relative reception overhead `epsilon`, i.e.,
`epsilon=0.2` is a 20% reception overhead.
"""
function ltfailure_lower(k::Int, epsilon::Number, Omega::Distribution{Univariate, Discrete}) :: Float64
    r = zero(Float64)
    tiny = realmin(Float64)
    for i in 1:div(k, 2)
        v = zero(r)
        v += exp(logbinomial(k, 2i-1) + ltfailure_lower_inner(2i-1, k, epsilon, Omega))
        v -= exp(logbinomial(k, 2i) + ltfailure_lower_inner(2i, k, epsilon, Omega))

        # the terms become smaller. stop after reaching the smallest normally
        # represented float.
        if v < tiny
            break
        end        
        r += v
    end
    return r
end
