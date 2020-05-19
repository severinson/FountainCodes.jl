using FountainCodes, Test, LinearAlgebra, Distributions

function init(::Type{GF256}, K; M=K-1, δ=1e-6)
    dd = FountainCodes.Soliton(K, M, δ)
    lt = FountainCodes.LTQ(K, dd)
    d = FountainCodes.Decoder(lt)
    src = [Vector{GF256}([i % 256]) for i in 1:K]
    return lt, d, src
end

function init(::Type{Float64}, K; M=K-1, δ=1e-6)
    dd = FountainCodes.Soliton(K, M, δ)
    lt = FountainCodes.LTQ{Float64}(K, dd)
    d = FountainCodes.Decoder(lt)
    src = randn(K)
    return lt, d, src
end

function init(::Type{Vector{Float64}}, K; M=K-1, δ=1e-6)
    dd = FountainCodes.Soliton(K, M, δ)
    lt = FountainCodes.LTQ{Float64}(K, dd)
    d = FountainCodes.Decoder(lt)
    src = [Vector{Float64}([i]) for i in 1:K]
    return lt, d, src
end

compare(a, b) = a == b
compare(a::Float64, b::Float64) = isapprox(a, b, rtol=1e-3)
compare(a::Vector{Float64}, b::Vector{Float64}) = isapprox(a, b, rtol=1e-3)

function test_mvnormal(K=10, r=round(Int, K*1.3))
    lt, d, src = init(Float64, K)
    Vs_μ = [get_value(lt, X, src) for X in 1:r]
    Vs_Σ = diagm(0=>ones(r))
    Vs = CodedMvNormal(Vs_μ, Vs_Σ)
    dec = decode(lt, 1:r, Vs)
    for i in 1:K
        if !compare(dec[i], src[i])
            error("expected dec[$i] to be $(src[i]), but got $(dec[i])")
        end
    end
    return true
end
@test test_mvnormal()

"""test that get_constraint returns the same constraint at each call"""
function test_constraint(VT, K=10, r=K)
    lt, _, src = init(VT, K)
    constraints = [get_constraint(lt, X) for X in 1:r]
    for X in 1:r
        correct = constraints[X]
        constraint = get_constraint(lt, X)
        if constraint != correct
            error("expected the $X-th constraint to be $correct, but got $constraint")
        end
    end
    return true
end
@test test_constraint(GF256, 10)
@test test_constraint(GF256, 100)
@test test_constraint(GF256, 1000)

function test_encode(VT, K=10)
    lt, _, src = init(VT, K)
    covered = zeros(Int, K)
    for X in 1:K
        constraint = get_constraint(lt, X)
        covered[constraint.nzind] .+= 1
    end
    for (i, v) in enumerate(covered)
        if iszero(v)
            error("$i-th source symbol not covered")
        end
    end
    return true
end
@test test_encode(GF256, 10)
@test test_encode(GF256, 100)
@test test_encode(GF256, 1000)

"""test that LT and LTQ ltgenerate choose the same indices"""
function test_ltgenerate_1(K=1024, M=K-1, δ=1e-6)
    dd = FountainCodes.Soliton(K, M, δ)
    # p_lt, _, C_lt = init(1024)

    lt = LT(K, dd)
    ltq = LTQ(K, dd)
    C = [Vector{GF256}([i % 256]) for i in 1:lt.L]

    # p_ltq, _, C_ltq = init_gf256(1024)
    for i in 1:K+300
        s_lt = FountainCodes.ltgenerate(C, i, lt)
        s_ltq = FountainCodes.ltgenerate(C, i, ltq)
        if s_lt.neighbours != s_ltq.neighbours
            error("$(s_lt.neighbours) != $(s_ltq.neighbours) for ESI $i")
        end
    end
    return true
end
# @test test_ltgenerate_1()

"""test diagonalization"""
function test_diagonalize(VT, K=10, r=round(Int, K*1.3))
    lt, d, src = init(VT, K)
    Vs = [get_value(lt, X, src) for X in 1:r]
    for X in 1:r
        FountainCodes.add!(d, get_constraint(lt, X))
    end
    FountainCodes.diagonalize!(d, Vs)
    for i in 1:d.num_decoded
        rpi = d.rowperm[i]
        cpi = d.colperm[i]
        coef = d.sparse[rpi][cpi]
        correct = coef.*src[cpi]
        for ci in 1:K
            cpj = d.colperm[ci]
            coef = FountainCodes.getdense(d, rpi, cpj)
            if !iszero(coef)
                correct = correct .+ coef.*src[cpj]
            end
        end
        if !compare(Vs[rpi], correct)
            error("expected values[$rpi] to be $correct, but got $(d.values[rpi])")
        end
    end
    return true
end
@test test_diagonalize(GF256, 10)
@test test_diagonalize(GF256, 100)
@test test_diagonalize(GF256, 1000)
@test test_diagonalize(Float64, 10)
@test test_diagonalize(Float64, 100)
@test test_diagonalize(Float64, 1000)
@test test_diagonalize(Vector{Float64}, 10)
@test test_diagonalize(Vector{Float64}, 100)
@test test_diagonalize(Vector{Float64}, 1000)

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
            error("expected values[$rpi] to be $correct, but got $(d.values[rpi])")
        end
    end
    return true
end
@test test_solve_dense(GF256, 10)
@test test_solve_dense(GF256, 100)
@test test_solve_dense(GF256, 1000)
@test test_solve_dense(Float64, 10)
@test test_solve_dense(Float64, 100)
@test test_solve_dense(Float64, 1000)
@test test_solve_dense(Vector{Float64}, 10)
@test test_solve_dense(Vector{Float64}, 100)
@test test_solve_dense(Vector{Float64}, 1000)

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
@test test_backsolve(GF256, 10)
@test test_backsolve(GF256, 100)
@test test_backsolve(GF256, 1000)
@test test_backsolve(Float64, 10)
@test test_backsolve(Float64, 100)
@test test_backsolve(Float64, 1000)
@test test_backsolve(Vector{Float64}, 10)
@test test_backsolve(Vector{Float64}, 100)
@test test_backsolve(Vector{Float64}, 1000)

"""test that decoding succeeds"""
function test_decode(VT, K=10, r=round(Int, K*1.3))
    lt, _, src = init(VT, K)
    Random.seed!(K)
    Xs = sample(1:100r, r, replace=false)
    Vs = [get_value(lt, X, src) for X in Xs]
    dec = decode(lt, Xs, Vs)
    for i in 1:K
        if !compare(dec[i], src[i])
            error("expected dec[$i] to be $(src[i]), but got $(dec[i])")
        end
    end
    return true
end
@test test_decode(GF256, 10, 13)
@test test_decode(GF256, 100, 120)
@test test_decode(GF256, 1000, 1300)
@test test_decode(Float64, 10)
@test test_decode(Float64, 100)
@test test_decode(Float64, 1000, 1350)
@test test_decode(Vector{Float64}, 10)
@test test_decode(Vector{Float64}, 100)
@test test_decode(Vector{Float64}, 700)
