using RaptorCodes, Base.Test
function logfactorial_exact_test_1(low=1, high=20)
    for n in low:high
        r = RaptorCodes.logfactorial_exact(n)
        c = log(factorial(n))
        err = abs.(r-c)
        if err > 1e-10
            error("logfactorial($n) = $r != $c, i.e., an error of $err")
        end
    end
    return true
end
@test logfactorial_exact_test_1()

function logfactorial_exact_test_2(low=1, high=20)
    for n in low:high
        for k in 1:n
            r = RaptorCodes.logfactorial_exact(n, k)
            c = log(factorial(n)) - log(factorial(k))
            err = abs.(r-c)
            if err > 1e-10
                error("logfactorial($n,$k) = $r != $c, i.e., an error of $err")
            end
        end
    end
    return true
end
@test logfactorial_exact_test_2()

function logfactorial_approx_test_1(low=20, high=100)
    for i in low:high
        r = RaptorCodes.logfactorial_approx(i)
        c = RaptorCodes.logfactorial_exact(i)
        err = abs.(r-c)
        if err > 1e-10
            error("logfactorial_approx($i) = $r != $c, i.e., an error of $err")
        end
    end
    return true
end
@test logfactorial_approx_test_1()

function logfactorial_approx_test_2(low=20, high=100)
    for n in low:high
        for k in low:n
            r = RaptorCodes.logfactorial_approx(n, k)
            c = RaptorCodes.logfactorial_exact(n, k)
            err = abs.(r-c)
            if err > 1e-6
                error("logfactorial_approx($n,$k) = $r != $c, i.e., an error of $err")
            end
        end
    end
    return true
end
@test logfactorial_approx_test_2()

function logfactorial_approx_test_3(n=134000, k=div(134000,2))
    r = RaptorCodes.logfactorial_approx(n, k)
    if isnan(r) || isinf(r)
        error("logfactorial_approx($n,$k) = $r")
    end
    return true
end
@test logfactorial_approx_test_3()

function logfactorial_approx_test_3(n=133999, k=2927)
    r = RaptorCodes.logfactorial_approx(n, k)
    if isnan(r) || isinf(r)
        error("logfactorial_approx($n,$k) = $r")
    end
    return true
end
@test logfactorial_approx_test_3()

function logfactorial_test_1(low=1, high=200)
    for n in low:high
        for k in low:n
            r = RaptorCodes.logfactorial(n, k)
            c = RaptorCodes.logfactorial_exact(n, k)
            err = abs.(r-c)
            if err > 1e-3
                error("logfactorial($n,$k) = $r != $c, i.e., an error of $err")
            end
        end
    end
    return true
end
@test logfactorial_test_1()

function logbinomial_test_1(low=1, high=20)
    for n in low:high
        for k in 1:n
            r = RaptorCodes.logbinomial(n, k)
            c = log(binomial(n, k))
            err = abs.(r-c)
            if err > 1e-6
                error("logbinomial($n,$k) = $r != $c, i.e., an error of $err")
            end
        end
    end
    return true
end
@test logbinomial_test_1()

function logbinomial_test_2(n=133999, k=2927)
    r = RaptorCodes.logbinomial(n, k)
    if isnan(r) || isinf(r)
        error("logbinomial($n,$k) = $r")
    end
    return true
end
@test logbinomial_test_2()

# k = 10000, M = 142, delta = 0.0317

function ltfailure_lower_test_1(low=1, high=20)
    for k in low:high
        M = Int(round(2/3*k))
        delta = 0.01
        Omega = Soliton(k, M, delta)
        for epsilon in linspace(0, 0.5, 10)
            r = ltfailure_lower(k, epsilon, Omega)
            c = RaptorCodes.ltfailure_lower_reference(k, epsilon, Omega)
            err = abs.(r-c)
            if err > 1e-3
                error("ltfailure_lower($k, $epsilon) = $r != $c, i.e., an error of $err")
            end
        end
    end
    return true
end
@test ltfailure_lower_test_1()

function ltfailure_lower_test_2()
    k = 10000
    M = 139
    delta = 0.0317
    Omega = Soliton(k, M, delta)
    println("mean=$(mean(Omega))")
    for (epsilon, c) in zip([0, 0.05, 0.1, 0.15], [1e-2, 1e-2/6, 1e-2/9, 1e-3])
        r = ltfailure_lower(k, epsilon, Omega)
        println("(epsilon, c)=($epsilon, $c), Omega.c=$(Omega.c), r=$r")
        # err = abs.(r-c)
        # if err > 1e-3
        #     println("ltfailure_lower($k, $epsilon) = $r")
        #     # println("ltfailure_lower($k, $epsilon) = $r != $c, i.e., an error of $err")
        # end
    end
    return false
end
# @test ltfailure_lower_test_2()

function ltfailure_lower_test_3()
    K = 4000
    M = 3998
    delta = 0.9999999701976676
    Omega = Soliton(K, M, delta)
    failure = [ltfailure_lower(K, epsilon, Omega) for epsilon in linspace(0.25, 0.4, 10)]
    println(failure)
    return false
end
# @test ltfailure_lower_test_3()
