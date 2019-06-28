using FountainCodes, Test, Distributions

function init_gf256(k=10)
    dd = FountainCodes.Soliton(k, round(Int, k*2/3), 0.01)
    p = FountainCodes.LTQ(k, dd)
    d = FountainCodes.Decoder(p)
    C = [Vector{GF256}([i % 256]) for i in 1:p.L]
    return p, d, C
end

function init_float64(k=10)
    dd = FountainCodes.Soliton(k, round(Int, k*2/3), 0.01)
    p = FountainCodes.LTQ{Float64}(k, dd)
    d = FountainCodes.Decoder(p)
    C = [Vector{Float64}([i]) for i in 1:p.L]
    return p, d, C
end

function init_float64_scalar(k=10)
    dd = FountainCodes.Soliton(k, round(Int, k*2/3), 0.01)
    p = FountainCodes.LTQ{Float64}(k, dd)
    d = FountainCodes.Decoder{Float64}(p)
    C = randn(p.L)
    return p, d, C
end

"test encoder using GF256 coefficients"
function test_encode_gf256()
    p, _, C = init_gf256()
    for i in 1:length(C)
        if !isassigned(C, i)
            error("intermediate symbol at index $i not assigned.")
        end
    end
    s = FountainCodes.ltgenerate(C, 1, p)
    return true
end
@test test_encode_gf256()

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

function test_subtract_gf256_1()
    p, d, C = init_gf256(10)
    s = FountainCodes.ltgenerate(C, 1, p)
    FountainCodes.add!(d, s)
    FountainCodes.add!(d, s)
    FountainCodes.expand_dense!(d)
    FountainCodes.subtract!(d, 1, 2, GF256(1))
    if !iszero(d.values[2])
        error("values[2]=$(d.values[2]) should be zero")
    end
    return true
end
@test test_subtract_gf256_1()

function test_subtract_gf256_2()
    p, d, C = init_gf256(10)
    a = GF256(3)
    FountainCodes.add!(d, QSymbol(1, [a], [1, 2], Vector{GF256}([1, 2])))
    FountainCodes.add!(d, QSymbol(1, [GF256(2)*a], [1, 2], Vector{GF256}([2, 4])))
    FountainCodes.expand_dense!(d)
    FountainCodes.subtract!(d, 1, 2, GF256(2))
    if !iszero(d.values[2])
        error("values[2]=$(d.values[2]) should be zero")
    end
    return true
end
@test test_subtract_gf256_2()

function test_subtract_float64_1()
    p, d, C = init_float64(10)
    s = FountainCodes.ltgenerate(C, 1, p)
    FountainCodes.add!(d, s)
    FountainCodes.add!(d, s)
    FountainCodes.expand_dense!(d)
    FountainCodes.subtract!(d, 1, 2, Float64(1))
    if !iszero(d.values[2])
        error("values[2]=$(d.values[2]) should be zero")
    end
    return true
end
@test test_subtract_float64_1()

function test_zerodiag_gf256()
    p, d, C = init_gf256(10)
    FountainCodes.add!(d, QSymbol(1, [GF256(2)], [1], Vector{GF256}([2])))
    FountainCodes.add!(d, QSymbol(1, [GF256(3)], [1], Vector{GF256}([3])))
    d.num_decoded += 1
    FountainCodes.expand_dense!(d)
    FountainCodes.zerodiag!(d, 2)
    if !iszero(d.values[2])
        error("values[2]=$(d.values[2]) should be zero")
    end
    return true
end
@test test_zerodiag_gf256()

function test_inactivate_gf256()
    p, d, C = init_gf256(10)
    FountainCodes.add!(d, QSymbol(1, [GF256(2)], [1], Vector{GF256}([2])))
    FountainCodes.add!(d, QSymbol(1, [GF256(3)], [1], Vector{GF256}([3])))
    d.num_decoded += 1 # simulate diagonalize! running for 1 iteration
    FountainCodes.mark_inactive!(d, 1)
    FountainCodes.setinactive!(d, 1)
    FountainCodes.setinactive!(d, 2)
    c = FountainCodes.getdense(d, 1, 1)
    if c != GF256(2)
        error("inactivated element of row 1 is $c but should be $(GF256(2))")
    end
    c = FountainCodes.getdense(d, 2, 1)
    if c != GF256(3)
        error("inactivated element of row 2 is $c but should be $(GF256(3))")
    end
    FountainCodes.subtract!(d, 1, 2, GF256(3)/GF256(2))
    c = FountainCodes.getdense(d, 2, 1)
    if !iszero(c)
        error("inactivated element of row 2 is $c but should be zero")
    end
    return true
end
@test test_inactivate_gf256()

function test_diagonalize_gf256_1()
    p, d, C = init_gf256(10)
    for i in 1:15
        s = FountainCodes.ltgenerate(C, i, p)
        FountainCodes.add!(d, s)
    end
    FountainCodes.diagonalize!(d)
    for i in 1:d.num_decoded
        rpi = d.rowperm[i]
        cpi = d.colperm[i]
        coef = d.sparse[rpi][cpi]
        correct = coef.*C[cpi]
        for ci in 1:d.p.L
            cpj = d.colperm[ci]
            coef = FountainCodes.getdense(d, rpi, cpj)
            if !iszero(coef)
                correct = correct + coef .* C[cpj]
            end
        end
        if d.values[rpi] != correct
            error("diagonalization failed. values[$rpi] is $(d.values[rpi]) but should be $correct")
        end
    end
    return true
end
@test test_diagonalize_gf256_1()

function test_diagonalize_gf256_2()
    p, d, C = init_gf256(20)
    for i in 1:22
        s = FountainCodes.ltgenerate(C, i, p)
        FountainCodes.add!(d, s)
    end
    FountainCodes.diagonalize!(d)
    for i in 1:d.num_decoded
        rpi = d.rowperm[i]
        cpi = d.colperm[i]
        coef = d.sparse[rpi][cpi]
        correct = coef.*C[cpi]
        for ci in 1:d.p.L
            cpj = d.colperm[ci]
            coef = FountainCodes.getdense(d, rpi, cpj)
            if !iszero(coef)
                correct = correct + coef .* C[cpj]
            end
        end
        if d.values[rpi] != correct
            error("diagonalization failed. values[$rpi] is $(d.values[rpi]) but should be $correct")
        end
    end
    return true
end
@test test_diagonalize_gf256_2()

function test_diagonalize_float64_1()
    p, d, C = init_float64(100)
    for i in 1:110
        s = FountainCodes.ltgenerate(C, i, p)
        FountainCodes.add!(d, s)
    end
    FountainCodes.diagonalize!(d)
    for i in 1:d.num_decoded
        rpi = d.rowperm[i]
        cpi = d.colperm[i]
        coef = d.sparse[rpi][cpi]
        correct = coef.*C[cpi]
        for ci in 1:d.p.L
            cpj = d.colperm[ci]
            coef = FountainCodes.getdense(d, rpi, cpj)
            if !iszero(coef)
                correct = correct + coef .* C[cpj]
            end
        end
        if !isapprox(d.values[rpi], correct, rtol=1e-3)
            err = abs(d.values[rpi] - correct)
            error("diagonalization failed. values[$rpi]=$(d.values[rpi]) but should be $correct. error is $err")
        end
    end
    return true
end
@test test_diagonalize_float64_1()

function test_ge_gf256_1()
    p, d, C = init_gf256()
    for i in 1:15
        s = FountainCodes.ltgenerate(C, i, p)
        FountainCodes.add!(d, s)
    end
    FountainCodes.diagonalize!(d)
    FountainCodes.solve_dense!(d)
    for i in d.p.L-d.num_inactivated+1:d.p.L
        rpi = d.rowperm[i]
        cpi = d.colperm[i]
        coef = FountainCodes.getdense(d, rpi, cpi)
        correct = coef.*C[cpi]
        if d.values[rpi] != correct
            error("GE failed. values[$rpi] is $(d.values[rpi]) but should be $correct")
        end
        if FountainCodes.countnz(d.dense, rpi) != 1
            error("GE failed. row[$rpi]=$(getcolumn(d.dense,rpi)) does not sum to 1.")
        end
    end
    return true
end
@test test_ge_gf256_1()

function test_ge_gf256_2()
    p, d, C = init_gf256(20)
    for i in 1:25
        s = FountainCodes.ltgenerate(C, i, p)
        FountainCodes.add!(d, s)
    end
    FountainCodes.diagonalize!(d)
    FountainCodes.solve_dense!(d)
    for i in d.p.L-d.num_inactivated+1:d.p.L
        rpi = d.rowperm[i]
        cpi = d.colperm[i]
        coef = FountainCodes.getdense(d, rpi, cpi)
        correct = coef.*C[cpi]
        if d.values[rpi] != correct
            error("GE failed. values[$rpi] is $(d.values[rpi]) but should be $correct")
        end
        if FountainCodes.countnz(d.dense, rpi) != 1
            error("GE failed. row[$rpi]=$(getcolumn(d.dense,rpi)) does not sum to 1.")
        end
    end
    return true
end
@test test_ge_gf256_2()

function test_backsolve_gf256()
    p, d, C = init_gf256(20)
    for i in 1:25
        s = FountainCodes.ltgenerate(C, i, p)
        FountainCodes.add!(d, s)
    end
    FountainCodes.diagonalize!(d)
    FountainCodes.solve_dense!(d)
    FountainCodes.backsolve!(d)
    for ri in 1:d.p.L-d.num_inactivated
        rpi = d.rowperm[ri]
        for ci in d.p.L-d.num_inactivated+1:d.p.L
            cpi = d.colperm[ci]
            coef = FountainCodes.getdense(d, rpi, cpi)
            if !iszero(coef)
                error("backsolve failed. row $ri column $ci is non-zero.")
            end
        end
    end
    return true
end
@test test_backsolve_gf256()

function test_decoder_gf256_1()
    p, d, C = init_gf256()
    for i in 1:15
        s = FountainCodes.ltgenerate(C, i, p)
        FountainCodes.add!(d, s)
    end
    output = FountainCodes.decode!(d)
    for i in 1:p.K
        if output[i] != C[i]
            error("decoding failure. source[$i] is $(output[i]). should be $(C[i]).")
        end
    end
    return true
end
@test test_decoder_gf256_1()

function test_decoder_gf256_2()
    p, d, C = init_gf256(1000)
    for i in 1:1400
        s = FountainCodes.ltgenerate(C, i, p)
        FountainCodes.add!(d, s)
    end
    output = FountainCodes.decode!(d)
    for i in 1:p.K
        if output[i] != C[i]
            error("decoding failure. source[$i] is $(output[i]). should be $(C[i]).")
        end
    end
    return true
end
@test test_decoder_gf256_2()

function test_decoder_gf256_3()
    K = 4000
    mode = 3998
    delta = 0.9999999701976676
    dd = FountainCodes.Soliton(K, mode, delta)
    p = LTQ(K, dd)
    d = FountainCodes.Decoder(p)
    C = [Vector{GF256}([i % 256]) for i in 1:p.L]
    for i in 1:6000
        s = FountainCodes.ltgenerate(C, i, p)
        FountainCodes.add!(d, s)
    end
    output = FountainCodes.decode!(d)
    for i in 1:p.K
        if output[i] != C[i]
            error("decoding failure. source[$i] is $(output[i]). should be $(C[i]).")
        end
    end
    return true
end
@test test_decoder_gf256_3()

function test_decoder_gf256_4()
    p, d, C = init_gf256(1024)
    for i in 1:1400
        s = FountainCodes.ltgenerate(C, i, p)
        FountainCodes.add!(d, s)
    end
    output = FountainCodes.decode!(d)
    for i in 1:p.K
        if output[i] != C[i]
            error("decoding failure. source[$i] is $(output[i]). should be $(C[i]).")
        end
    end
    return true
end
@test test_decoder_gf256_4()

function test_decoder_gf256_5()
    p, d, C = init_gf256(100)
    for i in 1:120
        s = FountainCodes.ltgenerate(C, i, p)
        FountainCodes.add!(d, s)
    end
    output = FountainCodes.decode!(d)
    for i in 1:p.K
        if output[i] != C[i]
            error("decoding failure. source[$i] is $(output[i]). should be $(C[i]).")
        end
    end
    return true
end
@test test_decoder_gf256_5()

"test decoding a dense q-ary LT code"
function test_dense_GF256()
    K = 100
    dd = DiscreteUniform(K/2, K)
    p = FountainCodes.LTQ(K, dd)
    d = FountainCodes.Decoder(p)
    C = [Vector{GF256}([i % 256]) for i in 1:p.L]
    for i in 1:K+10
        s = FountainCodes.ltgenerate(C, i, p)
        FountainCodes.add!(d, s)
    end
    output = FountainCodes.decode!(d)
    for i in 1:p.K
        if output[i] != C[i]
            error("decoding failure. source[$i] is $(output[i]). should be $(C[i]).")
        end
    end
    return true
end
@test test_dense_GF256()

"test decoding a dense LT code over the reals"
function test_dense_float64()
    K = 100
    dd = DiscreteUniform(K/2, K)
    p = FountainCodes.LTQ{Float64}(K, dd)
    d = FountainCodes.Decoder(p)
    C = [Vector{Float64}([i]) for i in 1:p.L]
    for i in 1:K+10
        s = FountainCodes.ltgenerate(C, i, p)
        FountainCodes.add!(d, s)
    end
    output = FountainCodes.decode!(d)
    for i in 1:p.K
        if !isapprox(output[i], C[i], rtol=1e-3)
            err = abs.(output[i] - C[i])
            error("decoding failure. source[$i] is $(output[i]). should be $(C[i]). error is $err.")
        end
    end
    return true
end
@test test_dense_float64()

function test_decoder_float64_scalar_1()
    p, d, C = init_float64_scalar()
    for i in 1:15
        s = FountainCodes.ltgenerate(C, i, p)
        FountainCodes.add!(d, s)
    end
    output = FountainCodes.decode!(d)
    for i in 1:p.K
        if !isapprox(output[i], C[i], rtol=1e-3)
            err = abs.(output[i] - C[i])
            error("decoding failure. source[$i] is $(output[i]). should be $(C[i]). error is $err.")
        end
    end
    return true
end
@test test_decoder_float64_scalar_1()

function test_decoder_float64_scalar_2()
    p, d, C = init_float64_scalar(100)
    for i in 1:130
        s = FountainCodes.ltgenerate(C, i, p)
        FountainCodes.add!(d, s)
    end
    output = FountainCodes.decode!(d)
    for i in 1:p.K
        if !isapprox(output[i], C[i], rtol=1e-3)
            err = abs.(output[i] - C[i])
            error("decoding failure. source[$i] is $(output[i]). should be $(C[i]). error is $err.")
        end
    end
    return true
end
@test test_decoder_float64_scalar_2()

function test_decoder_float64_1()
    p, d, C = init_float64()
    for i in 1:15
        s = FountainCodes.ltgenerate(C, i, p)
        FountainCodes.add!(d, s)
    end
    output = FountainCodes.decode!(d)
    for i in 1:p.K
        if !isapprox(output[i], C[i], rtol=1e-3)
            err = abs.(output[i] - C[i])
            error("decoding failure. source[$i] is $(output[i]). should be $(C[i]). error is $err.")
        end
    end
    return true
end
@test test_decoder_float64_1()

function test_decoder_float64_2()
    p, d, C = init_float64(100)
    for i in 1:130
        s = FountainCodes.ltgenerate(C, i, p)
        FountainCodes.add!(d, s)
    end
    output = FountainCodes.decode!(d)
    for i in 1:p.K
        if !isapprox(output[i], C[i], rtol=1e-3)
            err = abs.(output[i] - C[i])
            error("decoding failure. source[$i] is $(output[i]). should be $(C[i]). error is $err.")
        end
    end
    return true
end
@test test_decoder_float64_2()
