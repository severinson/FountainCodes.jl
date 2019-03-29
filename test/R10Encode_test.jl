using FountainCodes, Test

function init(k=10)
    p = FountainCodes.R10(k)
    C = [Vector{GF256}([i % 256]) for i in 1:p.L]
    N = [Dict{Int,Bool}() for _ in 1:p.L]
    precode!(C, p, N)
    return p, C, N
end

"make sure the encoder runs at all"
function test_encode_1()
    p, C, _ = init()
    for i in 1:length(C)
        if !isassigned(C, i)
            error("intermediate symbol at index $i not assigned.")
        end
    end
    s = FountainCodes.ltgenerate(C, 1, p)
    deg = FountainCodes.degree(s)
    if deg != 2
        error("LT degree is $deg bout should be 2")
    end
    return true
end
@test test_encode_1()

"test that the encoder produces a correct R10 constraint matrix"
function test_encode_2()
    p, C, N = init(10)
    ri = p.K+1
    correct = [1, 6, 7, 8]
    indices = sort!(collect(keys(N[ri])))
    if indices != correct
        error("R10 LDPC encoder failed. row $ri is $indices but should be $correct")
    end

    ri = p.K+7
    correct = [5, 6, 7, 10]
    indices = sort!(collect(keys(N[ri])))
    if indices != correct
        error("R10 LDPC encoder failed. row $ri is $indices but should be $correct")
    end

    ri = p.K+8
    correct = [1, 2, 4, 5, 8, 10, 11, 15]
    indices = sort!(collect(keys(N[ri])))
    if indices != correct
        error("R10 HDPC encoder failed. row $ri is $indices but should be $correct")
    end

    ri = p.K+11
    correct = [2, 3, 4, 5, 6, 7, 14, 15, 16, 17]
    indices = sort!(collect(keys(N[ri])))
    if indices != correct
        error("R10 HDPC encoder failed. row $ri is $indices but should be $correct")
    end
    return true

    ri = p.K+13
    correct = [11, 12, 13, 14, 15, 16, 17]
    indices = sort!(collect(keys(N[ri])))
    if indices != correct
        error("R10 HDPC encoder failed. row $ri is $indices but should be $correct")
    end
    return true
end
@test test_encode_2()
