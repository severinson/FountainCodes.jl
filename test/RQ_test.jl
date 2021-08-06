using FountainCodes, Test, SparseArrays

function init(K)
    error("Deprecated")
    p = FountainCodes.RQ(K)
    d = FountainCodes.Decoder(p)
    C = [Vector{GF256}([i % 256]) for i in 1:p.L]
    precode!(C, p)
    return p, d, C
end

# test RQ parameter choice
@test FountainCodes.RQ_parameters(27) == (30, 566, 11, 10, 41)
@test FountainCodes.RQ_parameters(30) == (30, 566, 11, 10, 41)
@test FountainCodes.RQ_parameters(90) == (91, 66, 17, 10, 103)

# test degree distribution
@test FountainCodes.RQ_deg(0) == 1
@test FountainCodes.RQ_deg(5243) == 2
@test FountainCodes.RQ_deg(1048576-1) == 30

"""test RaptorQ tuple function"""
function test_RQ_tuple()
    c = RQ(10)
    for X in 1:10
        d, a, b, d1, a1, b1 = FountainCodes.RQ_tuple(X, c)
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

"""test LDPC constraints against known correct constraints"""
function test_ldpc_constraints()
    c = RQ(10)
    A = FountainCodes.ldpc_constraint_matrix(c)

    ri = 1
    Is = [1, 6, 7, 8, 11, 18, 19]
    Vs = ones(GF256, 7)
    correct = sparsevec(Is, Vs, 27)
    ans = A[:, ri]
    @test ans == correct

    ri = 2
    Is = [1, 2, 7, 9, 12, 19, 20]
    Vs = ones(GF256, 7)
    correct = sparsevec(Is, Vs, 27)    
    ans = A[:, ri]
    @test ans == correct

    ri = 3
    Is = [1, 2, 3, 8, 10, 13, 20, 21]
    Vs = ones(GF256, 8)
    correct = sparsevec(Is, Vs, 27)        
    ans = A[:, ri]
    @test ans == correct        

    ri = 7
    Is = [5, 6, 7, 10, 17, 24, 25]
    Vs = ones(GF256, 7)
    correct = sparsevec(Is, Vs, 27)        
    ans = A[:, ri]
    @test ans == correct
    return
end
test_ldpc_constraints()

"""test HDPC constraints against known correct constraints"""
function test_hdpc_constraints()    
    c = RQ(10)
    A = FountainCodes.hdpc_constraint_matrix(c)

    i = 1
    Is = append!(collect(1:17), 18)    
    Vs = GF256.([
        0xfa, 0xf3, 0xf7, 0xf5, 0xf4, 0xf4, 0xf4, 0x7a, 0x3d,
        0x90, 0x48, 0x24, 0x12, 0x09, 0x04, 0x02, 0x01, 0x01,
    ])
    correct = sparsevec(Is, Vs, 27)    
    ans = A[:, i]
    @test ans == correct


    i = 10
    Is = append!(collect(1:17), 27)
    Vs = GF256.([
        0xeb, 0x75, 0x3a, 0x1d, 0x80, 0x40, 0x20, 0x10, 0x08,
        0x04, 0x02, 0x01, 0x8e, 0xc9, 0xea, 0x75, 0x3a ,0x01,
    ])
    correct = sparsevec(Is, Vs, 27)    
    ans = A[:, i]
    @test ans == correct
    return
end
test_hdpc_constraints()

"""Test that decoding with ESIs 0 through Kp-1 succeeds"""
function test_precode(K::Integer)
    rng = MersenneTwister(123)
    code = RQ(K)
    K = code.Kp # padding to the closest supported number of source symbols
    src = rand_nonzero(rng, GF256, K)

    # compute the intermediate symbols
    A = constraint_matrix(code, 0:K-1)
    pre = zeros(GF256, size(A, 1))
    pre[end-K+1:end] .= src
    pre = decode(A, pre)

    # test that the first K coded symbols are equal to the src symbols, i.e., that the code is 
    # systematic
    G = generator_matrix(code, 0:K-1)
    out = G'*pre
    @test out[1:K] == src
end
for K in 4:100
    test_precode(K)
end
test_precode(7000)
test_precode(8192)