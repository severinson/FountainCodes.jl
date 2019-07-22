using FountainCodes, BenchmarkTools

function init_lt(K; M=K-1, δ=1e-6)
    dd = FountainCodes.Soliton(K, M, δ)
    lt = FountainCodes.LT(K, dd)
    src = [Vector{GF256}([i % 256]) for i in 1:lt.L]
    return lt, src
end

function benchmark_lt(K; r=round(Int, K*1.3), M=K-1, δ=1e-9)
    println("LT Benchmark: K=$K, r=$r, M=$M, δ=$δ")
    lt, src = init_lt(K, M=M, δ=δ)
    Vs = [get_value(lt, X, src) for X in 1:r]
    decode(lt, 1:r, Vs) # warm-up the jit
    rv = @benchmark decode($lt, $(1:r), $Vs)
    show(stdout, "text/plain", rv)
    println("\n")
    return
end
benchmark_lt(1000, r=1300)
benchmark_lt(10000, r=13000)

function init_ltq(::Type{GF256}, K; M=K-1, δ=1e-6)
    dd = FountainCodes.Soliton(K, M, δ)
    lt = FountainCodes.LTQ(K, dd)
    src = [Vector{GF256}([i % 256]) for i in 1:lt.L]
    return lt, src
end

function benchmark_ltq(K; r=round(Int, K*1.3), M=K-1, δ=1e-9)
    println("LTQ Benchmark: K=$K, r=$r, M=$M, δ=$δ")
    lt, src = init_ltq(GF256, K, M=M, δ=δ)
    Vs = [get_value(lt, X, src) for X in 1:r]
    decode(lt, 1:r, Vs) # warm-up the jit
    rv = @benchmark decode($lt, $(1:r), $Vs)
    show(stdout, "text/plain", rv)
    println("\n")
    return
end
benchmark_ltq(1000, r=1300)
benchmark_ltq(10000, r=13000)
