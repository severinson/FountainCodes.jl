using FountainCodes, Test, Distributions, LinearAlgebra, SparseArrays

"""R10 code with scalar values"""
function init(::Type{GF256}, K)
    r10 = R10(K)
    d = Decoder{Bool}(K)
    src = [GF256(i % 256) for i in 1:K]
    inter = precode(src, r10)
    return r10, d, src, inter
end

"""R10 code with vector values"""
function init(::Type{Vector{GF256}}, K)
    r10 = R10(K)
    d = Decoder{Bool}(K)
    src = [Vector{GF256}([i % 256]) for i in 1:K]
    inter = precode(src, r10)
    return r10, d, src, inter
end

compare(a, b) = a == b

"""test the R10 degree distribution."""
function test_r10_degree()

    # test the values at which the degree increments
    f = [0, 10241, 491582, 712794, 831695, 948446, 1032189]
    correct = [1, 2, 3, 4, 10, 11, 40]
    for (v, d) in zip(f, correct)
        ans = FountainCodes.r10_degree(v)
        if ans != d error("expected r10_degree($v) to be $d, but got $ans.") end
    end

    # test the values one before the degree increments
    f = [10241, 491582, 712794, 831695, 948446, 1032189, 1048576]
    correct = [1, 2, 3, 4, 10, 11, 40]
    for (v, d) in zip(f, correct)
        ans = FountainCodes.r10_degree(v-1)
        if ans != d error("expected r10_degree($(v-1)) to be $d, but got $ans.") end
    end
    return true
end
@test test_r10_degree()

"""test that the constraint matrix for the first K symbols is of full rank"""
function test_systematic_indices(K=10)
    code = R10(K)
    Xs = -(code.S+code.H):(code.K-1)
    constraints = [get_constraint(code, X) for X in Xs]
    M = Matrix(hcat(constraints...))'
    @assert size(M) == (code.L, code.L)
    r = rank(M)
    if r != code.L
        error("expected constraint matrix to have rank $(code.L), but got $r.")
    end
    return true
end
for K in 4:100 @test test_systematic_indices(K) end
@test test_systematic_indices(1000)
@test test_systematic_indices(2500)
# @test test_systematic_indices(7000)

"""test the precode constraints"""
function test_precode_constraints()
    code = R10(10)

    # LDPC constraints
    X = -1
    correct = sparsevec([1, 6, 7, 8, 11], ones(Bool, 5), dimension(code))
    ans = get_constraint(code, X)
    if correct != ans
        error("expected 1st LDPC constraint to be\n$correct,\nbut got\n$ans")
    end

    X = -7
    correct = sparsevec([5, 6, 7, 10, 17], ones(Bool, 5), dimension(code))
    ans = get_constraint(code, X)
    if correct != ans
        error("expected 7th LDPC constraint to be\n$correct,\nbut got\n$ans")
    end

    # HDPC constraints
    X = -8
    Is = [1, 2, 4, 5, 8, 10, 11, 15, 18]
    Vs = ones(Bool, length(Is))
    correct = sparsevec(Is, Vs, dimension(code))
    ans = get_constraint(code, X)
    if correct != ans
        error("expected 1st HDPC constraint to be\n$correct,\nbut got\n$ans")
    end

    X = -11
    Is = [2, 3, 4, 5, 6, 7, 14, 15, 16, 17, 21]
    Vs = ones(Bool, length(Is))
    correct = sparsevec(Is, Vs, dimension(code))
    ans = get_constraint(code, X)
    if correct != ans
        error("expected 4th HDPC constraint to be\n$correct,\nbut got\n$ans")
    end

    X = -13
    Is = [11, 12, 13, 14, 15, 16, 17, 23]
    Vs = ones(Bool, length(Is))
    correct = sparsevec(Is, Vs, dimension(code))
    ans = get_constraint(code, X)
    if correct != ans
        error("expected 6th HDPC constraint to be\n$correct,\nbut got\n$ans")
    end
    return true
end
@test test_precode_constraints()

"""test that the precode values are computed correctly"""
function test_precode_values(K=10)
    code = R10(K)
    src = [GF256(i % 256) for i in 0:K-1]
    inter = precode(src, code)
    for X in 0:K-1
        correct = src[X+1]
        ans = get_value(code, X, inter)
        if correct != ans
            error("expected $X-th LT symbol to be $correct, but got $ans")
        end
    end
    return true
end
for K in 4:100 @test test_precode_values(K) end
@test test_precode_values(1000)
@test test_precode_values(2500)
# @test test_precode_values(7000)
# @test test_precode_values(8000)

"""test diagonalization"""
function test_diagonalize(VT, K=10, r=round(Int, K*1.3))
    code, d, src, inter = init(VT, K)
    Vs = [zero(inter[1]) for _ in 1:(code.S+code.H)] # parity symbol values
    append!(Vs, [get_value(code, X, inter) for X in 0:r-1])
    for X in -(code.S+code.H):r-1
        FountainCodes.add!(d, get_constraint(code, X))
    end
    FountainCodes.diagonalize!(d, Vs)
    for i in 1:d.num_decoded
        rpi = d.rowperm[i]
        cpi = d.colperm[i]
        coef = d.sparse[rpi][cpi]
        correct = coef.*inter[cpi]
        for ci in 1:K
            cpj = d.colperm[ci]
            coef = FountainCodes.getdense(d, rpi, cpj)
            if !iszero(coef)
                correct = correct .+ coef.*inter[cpj]
            end
        end
        if !compare(Vs[rpi], correct)
            error("expected Vs[$rpi] to be $correct, but got $(Vs[rpi])")
        end
    end
    return true
end
# @test test_diagonalize(GF256, 10)
# @test test_diagonalize(GF256, 100)
# @test test_diagonalize(GF256, 1000)
# @test test_diagonalize(Float64, 10)
# @test test_diagonalize(Float64, 100)
# @test test_diagonalize(Float64, 1000)
# @test test_diagonalize(Vector{Float64}, 10)
# @test test_diagonalize(Vector{Float64}, 100)
# @test test_diagonalize(Vector{Float64}, 1000)

"""test solving the dense subsystem u_lower"""
function test_solve_dense(VT, K=10, r=round(Int, K*1.3))
    lt, d, src = init(VT, K)
    Vs = [get_value(lt, X, src) for X in 1:r]
    for X in 1:r
        FountainCodes.add!(d, get_constraint(lt, X))
    end
    FountainCodes.diagonalize!(d, Vs)
    FountainCodes.solve_dense!(d, Vs)
    for i in K-d.num_inactivated+1:K
        rpi = d.rowperm[i]
        cpi = d.colperm[i]
        coef = FountainCodes.getdense(d, rpi, cpi)
        correct = coef.*src[cpi]
        if !compare(Vs[rpi], correct)
            error("expected Vs[$rpi] to be $correct, but got $(Vs[rpi])")
        end
    end
    return true
end
# @test test_solve_dense(GF256, 10)
# @test test_solve_dense(GF256, 100)
# @test test_solve_dense(GF256, 1000)
# @test test_solve_dense(Float64, 10)
# @test test_solve_dense(Float64, 100)
# @test test_solve_dense(Float64, 1000)
# @test test_solve_dense(Vector{Float64}, 10)
# @test test_solve_dense(Vector{Float64}, 100)
# @test test_solve_dense(Vector{Float64}, 1000)

"""test that backsolve zeroes out the dense matrix correctly"""
function test_backsolve(VT, K=10, r=round(Int, K*1.3))
    lt, d, src = init(VT, K)
    Vs = [get_value(lt, X, src) for X in 1:r]
    for X in 1:r
        FountainCodes.add!(d, get_constraint(lt, X))
    end
    FountainCodes.diagonalize!(d, Vs)
    FountainCodes.solve_dense!(d, Vs)
    FountainCodes.backsolve!(d, Vs)
    for ri in 1:K-d.num_inactivated
        rpi = d.rowperm[ri]
        for ci in K-d.num_inactivated+1:K
            cpi = d.colperm[ci]
            coef = FountainCodes.getdense(d, rpi, cpi)
            if !iszero(coef)
                error("expected dense[$ri, $ci] to be zero, but got $coef")
            end
        end
    end
    return true
end
# @test test_backsolve(GF256, 10)
# @test test_backsolve(GF256, 100)
# @test test_backsolve(GF256, 1000)
# @test test_backsolve(Float64, 10)
# @test test_backsolve(Float64, 100)
# @test test_backsolve(Float64, 1000)
# @test test_backsolve(Vector{Float64}, 10)
# @test test_backsolve(Vector{Float64}, 100)
# @test test_backsolve(Vector{Float64}, 1000)

"""test that decoding succeeds"""
function test_decode(VT, K=10, r=round(Int, K*1.3))
    code, d, src, inter = init(VT, K)
    Vs = [zero(inter[1]) for _ in 1:(code.S+code.H)] # parity symbol values
    append!(Vs, [get_value(code, X, inter) for X in 0:r-1]) # LT symbol values
    Xs = -(code.S+code.H):r-1 # ESIs incl. parity symbols
    dec_inter = decode(code, Xs, Vs)
    for i in 1:code.L
        if dec_inter[i] != inter[i]
            error("expected $(inter[i]), but got $(dec_inter[i])")
        end
    end
    return true
end
for K in 4:100 @test test_decode(GF256, K, K) end
@test test_decode(Vector{GF256}, 100, 100)
# @test test_decode(GF256, 8192, 8192)
