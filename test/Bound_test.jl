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

function logbinomial_test(low=1, high=20)
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
@test logbinomial_test()
