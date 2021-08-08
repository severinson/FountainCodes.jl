using FountainCodes, Test, SparseArrays

function test_r10_degree()

    # test the values at which the degree increments
    f = [0, 10241, 491582, 712794, 831695, 948446, 1032189]
    correct = [1, 2, 3, 4, 10, 11, 40]
    for (v, d) in zip(f, correct)
        ans = FountainCodes.r10_degree(v)
        if ans != d 
            error("expected r10_degree($v) to be $d, but got $ans.")
        end
    end

    # test the values one before the degree increments
    f = [10241, 491582, 712794, 831695, 948446, 1032189, 1048576]
    correct = [1, 2, 3, 4, 10, 11, 40]
    for (v, d) in zip(f, correct)
        ans = FountainCodes.r10_degree(v-1)
        if ans != d 
            error("expected r10_degree($(v-1)) to be $d, but got $ans.") 
        end
    end
    return true
end
@test test_r10_degree()

function test_ldpc_constraints()
    code = R10(10)
    A = FountainCodes.ldpc_constraint_matrix(code)

    # 1st LDPC constraint
    correct = sparsevec([1, 6, 7, 8, 11], ones(GF256, 5), dimension(code))
    ans = A[:, 1]
    @test ans == correct

    # 7th LDPC constraint
    correct = sparsevec([5, 6, 7, 10, 17], ones(GF256, 5), dimension(code))
    ans = A[:, 7]
    @test ans == correct
end
test_ldpc_constraints()

function test_hdpc_constraints()
    code = R10(10)
    A = FountainCodes.hdpc_constraint_matrix(code)    

    # 1st HDPC constraint
    Is = [1, 2, 4, 5, 8, 10, 11, 15, 18]
    Vs = ones(GF256, length(Is))
    correct = sparsevec(Is, Vs, dimension(code))
    ans = A[:, 1]
    @test ans == correct

    # 4th HDPC constraint
    Is = [2, 3, 4, 5, 6, 7, 14, 15, 16, 17, 21]
    Vs = ones(Bool, length(Is))
    correct = sparsevec(Is, Vs, dimension(code))
    ans = A[:, 4]
    @test ans == correct

    # 6th HDPC constraint
    Is = [11, 12, 13, 14, 15, 16, 17, 23]
    Vs = ones(Bool, length(Is))
    correct = sparsevec(Is, Vs, dimension(code))
    ans = A[:, 6]
    @test ans == correct
    return
end
test_hdpc_constraints()

"""Test that decoding with ESIs 0 through K-1 succeeds"""
function test_precode(K::Integer)
    rng = MersenneTwister(123)
    code = R10(K)
    src = rand_nonzero(rng, GF256, K)

    # compute the intermediate symbols
    A = constraint_matrix(code, 0:K-1)
    pre = zeros(GF256, size(A, 2))
    pre[end-K+1:end] .= src
    int = decode(A, pre)

    # test that the first K coded symbols are equal to the src symbols, i.e., that the code is 
    # systematic
    G = generator_matrix(code, 0:K-1)
    enc = G'*int
    @test enc[1:K] == src
end
for K in 4:100
    test_precode(K)
end
test_precode(7000)
test_precode(8192)