using FountainCodes, Test

"""return an LT code object and a vector of source symbols"""
function init(K; M=K-1, δ=1e-6)
    dd = FountainCodes.Soliton(K, M, δ)
    lt = FountainCodes.LT(K, dd)
    src = [Vector{GF256}([i % 256]) for i in 1:lt.L]
    return lt, src
end

"""test that get_constraint returns the same constraint at each call"""
function test_constraint(K=10, r=K)
    lt, src = init(K)
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
@test test_constraint(10)
@test test_constraint(100)
@test test_constraint(1000)

"""test that the encoder covers all source symbols"""
function test_encode(K=10)
    lt, src = init(K)
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
@test test_encode(10)
@test test_encode(100)
@test test_encode(1000)

"""test diagonalization"""
function test_diagonalize(K=10, r=round(Int, K*1.3))
    lt, src = init(K)
    Vs = [get_value(lt, X, src) for X in 1:r]
    d = Decoder(lt)
    for X in 1:r
        FountainCodes.add!(d, get_constraint(lt, X))
    end
    FountainCodes.diagonalize!(d, Vs)
    for i in 1:d.num_decoded
        rpi = d.rowperm[i]
        cpi = d.colperm[i]
        correct = src[cpi]
        for ci in 1:d.p.L
            cpj = d.colperm[ci]
            if !iszero(FountainCodes.getdense(d, rpi, cpj))
                correct = correct + src[cpj]
            end
        end
        if Vs[rpi] != correct
            error("expected values[$rpi] to be $correct, but got $(d.values[rpi])")
        end
    end
    return true
end
@test test_diagonalize(10)
@test test_diagonalize(100)
@test test_diagonalize(1000)

"""test solving the dense subsystem u_lower"""
function test_solve_dense(K=10, r=round(Int, K*1.3))
    lt, src = init(K)
    Vs = [get_value(lt, X, src) for X in 1:r]
    d = Decoder(lt)
    for X in 1:r
        FountainCodes.add!(d, get_constraint(lt, X))
    end
    FountainCodes.diagonalize!(d, Vs)
    FountainCodes.solve_dense!(d, Vs)
    for i in d.p.L-d.num_inactivated+1:d.p.L
        rpi = d.rowperm[i]
        cpi = d.colperm[i]
        correct = src[cpi]
        if Vs[rpi] != correct
            error("expected values[$rpi] to be $correct, but got $(d.values[rpi])")
        end
    end
    return true
end
@test test_solve_dense(10)
@test test_solve_dense(100)
@test test_solve_dense(1000)

"""test that decoding succeeds"""
function test_decode(K=10, r=round(Int, K*1.3))
    lt, src = init(K)
    Vs = [get_value(lt, X, src) for X in 1:r]
    dec = decode(lt, 1:r, Vs)
    for i in 1:lt.K
        if dec[i] != src[i]
            error("expected $(src[i]), but got $(dec[i])")
        end
    end
    return true
end
@test test_decode(10, 13)
@test test_decode(100, 120)
@test test_decode(1000, 1300)
