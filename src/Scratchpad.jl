# Development scratchpad

using Random
using StatsBase, Statistics, Distributions
using Primes, DataStructures, SparseArrays

# type system
abstract type CoefficientType end
mutable struct Binary <: CoefficientType end
mutable struct NonBinary <: CoefficientType end
abstract type Code{T<:CoefficientType} end
const BinaryCode = Code{Binary}
const NonBinaryCode = Code{NonBinary}
abstract type Selector end

"symbol value"
abstract type Value end

"coded symbol"
abstract type CodeSymbol end

include("Bound.jl")
include("Numinv.jl")
include("Soliton.jl")
include("GF256.jl")
include("CodedMvNormal.jl")
include("Symbols.jl")
include("Gray.jl")
include("QMatrix.jl")
include("DisjointSet.jl")
include("Decode.jl")
include("HeapSelect.jl")
include("R10Tables.jl")
include("R10.jl")
include("LT.jl")
include("RQTables.jl")
include("RQ.jl")
include("LDPC.jl")

# QMatrix Benchmarks

function xora!(c::BitArray, a::BitArray, b::BitArray, n::Int)::BitArray
    for _ in 1:n
        for i in 1:length(a.chunks)
            c.chunks[i] = xor(a.chunks[i], b.chunks[i])
        end
    end
    return c
end

function xorb!(c::BitArray, a::BitArray, b::BitArray, n::Int)::BitArray
    for _ in 1:n
        c .= xor.(a, b)
    end
    return c
end

"""compare xoring the chunks of bit arrays with xor."""
function benchmark_xor()
    n = 1000000
    m = 64*100
    a = trues(m)
    b = falses(m)
    c1 = falses(m)
    println("xora")
    @time c1 .= xora!(c1, a, b, n)
    c2 = falses(m)
    println("xorb")
    @time c2 .= xorb!(c2, a, b, n)
    @assert c1 == c2
    return
end

function subtractn(M, c1, c2, n)
    for _ in 1:n
        subtract!(M, c1, c2)
    end
end

"""benchmark the column subtract operation"""
function benchmark_subtract()
    m, n = 64*1000, 10000
    M = QMatrix{Float64}(m, n)
    M[:, 1] .= 0.1
    M[:, 2] .= true
    # subtract!(M, 2.0, 1, 2)
    @timev subtractn(M, 1, 2, 10000)
    return
end

# GF256 arithmetic benchmarks

function subtract_subeq!(a, b, c, n)
    for _ in 1:n
        subeq!(a, b, c)
    end
end

function subtract_dot!(a, b, c, n)
    for _ in 1:n
        a .-= c.*b
        println(a)
    end
end

function benchmark_subeq()
    c = rand(GF256)
    n, m = 1000, 100000
    # n, m = 100, 1000
    a = [rand(GF256) for _ in 1:m]
    b = [rand(GF256) for _ in 1:m]
    # println("subeq!")
    # subtract_subeq!(a, b, c, n)
    println("dot")
    subtract_dot!(a, b, c, 2)
end

# vdegree benchmarks

mutable struct D
    num_decoded::Int
    num_inactivated::Int
    num_symbols::Int
end

function vdegree_loop_inner(v, p, d)
    degree = 0
    i = d.num_decoded
    u = d.num_inactivated
    L = d.num_symbols
    # @inbounds return reduce((v1, v2)->v1 + i < v2 < L-u, v.nzind)
    # @views return sum(i .< p[v.nzind] .< L-u)
    # @inbounds for k in 1:length(v.nzind)
    @inbounds for j in v.nzind
        if i < p[j] < L-u
            # if i < p[v.nzind[k]] < L-u
            degree += 1
        end
    end
    return degree
end

function vdegree_loop(vs, p, d, n)
    total = 0
    for i in 1:n
        # v = vs[i % length(vs) + 1]
        # @inbounds total += vdegree_loop_inner(v, p, d)
        @inbounds total += vdegree_loop_inner(vs[i % length(vs) + 1], p, d)
        # @inbounds total += vdegree_loop_inner(vs[1], p, d)
    end
    return total
end

function vdegree_vector_inner(v, s)
    @views return sum(s[v.nzind])
end

function vdegree_vector(vs, s, n)
    total = 0
    for i in 1:n
        # v = vs[i % length(vs) + 1]
        # @inbounds total += vdegree_loop_inner(v, p, d)
        @inbounds total += vdegree_vector_inner(vs[i % length(vs) + 1], s)
        # @inbounds total += vdegree_loop_inner(vs[1], p, d)
    end
    return total
end

function benchmark_vdegree(n=10000, m=1000, nv=100, nsamples=1000)
    Is = sample(1:n, m, replace=false)
    Vs = ones(Int, m)
    vs = [sparsevec(Is, Vs, n) for _ in 1:nv]
    p = randperm(n)
    u = 100
    l = 100

    println("bitarray")
    s = trues(n)
    s[1:u] .= false
    s[end-l:end] .= false
    t = @elapsed @time vdegree_vector(vs, s, nsamples)
    t /= nv * nsamples
    return t

    # println("loop")
    # d = D(100, 100, n)
    # t = @elapsed @time vdegree_loop(vs, p, d, nsamples)
    # t /= nv * nsamples
    # return t
end

## callback benchmarks

function array_sample(i, j, v, n)
    for _ in 1:n v[i] -= v[j] end
end

function callback_sample(i, j, f, n)
    for _ in 1:n f(i, j) end
end

fun = () -> nothing

function global_callback_sample(i, j, n)
    for _ in 1:n fun(i, j) end
end

function benchmark_callback(n=10000, nsamples=100000)
    v = randn(n)

    println("array")
    @time array_sample(1, 10, v, nsamples)

    println("callback")
    @time callback_sample(2, 11, (i, j)->v[i]-=v[j], nsamples)

    println("callback sep.")
    f = (i, j)->v[i]-=v[j]
    @time callback_sample(2, 11, f, nsamples)

    println("global callback")
    global fun = (i, j)->v[i]-=v[j]
    @time global_callback_sample(2, 11, nsamples)

    return
end

using LinearAlgebra

## MvNormal ##
# function subtract!()

function mvnormal_main(K=10, r=round(Int, K*1.05), M=K-1, δ=1e-9)
    Random.seed!()
    # generate source data
    # c = rand(K, K) .+ diagm(0=>ones(K))
    c = diagm(0=>ones(K))
    c[1, 1] = 1
    c[1, 2] = -1
    # c[1:end, 1] .+= 1:K
    # c[1, 1:end] .+= 1:K
    # c = diagm(0=>ones(K))
    # c = Matrix(1.0.*I, K, K)
    # @assert rank(c) == K
    Σ_noise = c*Matrix(1.0.*I, K, K)*c'
    @assert isapprox(Σ_noise, Σ_noise')
    Noise = MvNormal(Σ_noise)
    # src = Float64.(1:K)
    src = randn(K)
    src_noise = src .+ rand(Noise)::Vector{Float64}

    # encode
    dd = Soliton(K, M, δ)
    lt = LTQ{Float64}(K, dd)
    C = Matrix(vcat([get_constraint(lt, X)' for X in 1:r]...))
    Σ_enc = C*Σ_noise*C'
    Vs = CodedMvNormal([get_value(lt, X, src_noise) for X in 1:r], Σ_enc)
    Random.seed!()
    @assert rank(C) == K
    @assert isapprox(C*src_noise, mean(Vs))

    # optimal estimate without coding (using PDMats transform)
    d_uncoded = MvNormal(src_noise, Σ_noise)
    M_X = PDMats.whiten(d_uncoded.Σ, Matrix(1.0.*I, K, K))
    M_src = PDMats.whiten(d_uncoded.Σ, src_noise)
    β = M_X \ M_src
    se = (β .- src).^2
    mse = mean(se)
    println("optimal:\t", mse)

    # optimal estimate without coding (using my transform)
    d = CodedMvNormal(src_noise, Σ_noise)
    M = whiten(d)
    β = M \ (M*src_noise)
    se = (β .- src).^2
    mse = mean(se)
    println("optimal (mine):\t", mse)

    # estimate without whitening
    β_ls = C \ mean(Vs)
    se_ls = (β_ls .- src).^2
    mse_ls = mean(se_ls)
    println("ls:\t\t", mse_ls)

    # estimate with whitening
    M = whiten(Vs)
    β_ls_wt = (M*C) \ (M*mean(Vs))
    se_ls_wt = (β_ls_wt .- src).^2
    mse_ls_wt = mean(se_ls_wt)
    println("ls (wt):\t", mse_ls_wt)

    # decoding
    dec = decode(lt, 1:r, deepcopy(mean(Vs)))
    se = (dec .- src).^2
    mse = mean(se)
    println("decoded:\t", mse)

    # decoding with whitening
    dec = decode(lt, 1:r, deepcopy(Vs))
    se = (dec .- src).^2
    mse = mean(se)
    println("decoded (wt):\t", mse)
end
