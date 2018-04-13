using RaptorCodes, Base.Test

function init(K=10)
    p = RaptorCodes.RQ(K)
    C = Vector{Vector{GF256}}(p.L)
    for i = 1:p.K
        C[i] = Vector{GF256}([i % 256])
    end
    precode!(C, p, N)
    return p, C, N
end

# test RQ parameter choice
@test RaptorCodes.RQ_parameters(27) == (30, 566, 11, 10, 41)
@test RaptorCodes.RQ_parameters(30) == (30, 566, 11, 10, 41)
@test RaptorCodes.RQ_parameters(90) == (91, 66, 17, 10, 103)

# test degree distribution
@test RaptorCodes.RQ_deg(0) == 1
@test RaptorCodes.RQ_deg(5243) == 2
@test RaptorCodes.RQ_deg(1048576-1) == 30

"""test tuple function"""
function test_RQ_tuple()
    c = RQ(10)
    for X in 1:10
        d, a, b, d1, a1, b1 = RaptorCodes.RQ_tuple(X, c)
        if !(0 < d)
            error("d must be positive")
        end
        if !(1 <= a <= c.W-1)
            error("a must be between 1 and W-1 inclusive")
        end
        if !(0 <= b <= c.W-1)
            error("b must be between 0 and W-1 inclusive")
        end
        if !(d1 in [2, 3])
            error("d1 must be either 2 or 3")
        end
        if !(1 <= a1 <= c.P1-1)
            error("a1 must be between 1 and P1-1 inclusive")
        end
        if !(0 <= b1 <= c.P1-1)
            error("b1 must be between 0 and P1-1 inclusive")
        end
    end
    return true
end
@test test_RQ_tuple()

"make sure the encoder runs at all"
function test_precode_relations()
    c = RQ(10)
    N = RaptorCodes.precode_relations(c)

    # test LDPC constraints
    ri = 1
    correct = [1, 6, 7, 8, 11, 18, 19]
    indices = N[ri][1]
    # indices = sort!(collect(keys(N[ri])))
    if indices != correct
        error("incorrect RQ LDPC constraint. row $ri is $indices but should be $correct")
    end

    ri = 7
    correct = [5, 6, 7, 10, 17, 24, 25]
    indices = N[ri][1]
    if indices != correct
        error("incorrect RQ LDPC constraint. row $ri is $indices but should be $correct")
    end

    # test HDPC constraints
    ri = 8
    correct = append!(collect(1:17), 18)
    indices = N[ri][1]
    if indices != correct
        error("RQ HDPC constraint error. row $ri indices are\n$indices,\nbut should be\n$correct")
    end
    correct = [
        0xfa, 0xf3, 0xf7, 0xf5, 0xf4, 0xf4, 0xf4, 0x7a, 0x3d,
        0x90, 0x48, 0x24, 0x12, 0x09, 0x04, 0x02, 0x01, 0x01,
    ]
    coefs = N[ri][2]
    if coefs != correct
        error("RQ HDPC constraint error. row $ri coefficients are\n$coefs,\n but should be\n$correct")
    end

    ri = 17
    correct = append!(collect(1:17), 27)
    indices = N[ri][1]
    if indices != correct
        error("RQ HDPC constraint error. row $ri indices are\n$indices,\nbut should be\n$correct")
    end
    correct = [
        0xeb, 0x75, 0x3a, 0x1d, 0x80, 0x40, 0x20, 0x10, 0x08,
        0x04, 0x02, 0x01, 0x8e, 0xc9, 0xea, 0x75, 0x3a ,0x01,
    ]
    coefs = N[ri][2]
    if coefs != correct
        error("RQ HDPC constraint error. row $ri coefficients are\n$coefs,\n but should be\n$correct")
    end
    return true
end
@test test_precode_relations()

function test_ltgenerate()
    c = RQ(10)
    C = zeros(GF256, c.L)
    for X in 1:100
        s = ltgenerate(C, X, c)
        for i in 2:length(s.neighbours)
            if s.neighbours[i] <= s.neighbours[i-1]
                error("indices for symbol $s are not strictly increasing")
            end
        end
        if s.neighbours[1] < 1 || s.neighbours[end] > c.L
            error("indices for symbol $s are not between 1 and L")
        end
    end
    return true
end
@test test_ltgenerate()

function test_precode()
    c = RQ(10)
    C = [Vector{GF256}([i]) for i in 1:c.L]
    println("C=$C")
    precode!(C, c)
    println("C=$C")
    for X in 1:c.K
        s = ltgenerate(C, X, c)
        correct = Vector{GF256}([X])
        if s.value != correct
            error("$X-th LT symbol value should be $correct but is $(s.value)")
        end
    end
    return true
end
@test test_precode()