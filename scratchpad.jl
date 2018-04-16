using RaptorCodes

function init(k=2000)
    p = RaptorCodes.R10Parameters(k)
    d = RaptorCodes.Decoder(p)
    C = Vector{Vector{F256}}(p.L)
    for i = 1:p.K
        C[i] = Vector{F256}([i % 256])
    end
    RaptorCodes.r10_ldpc_encode!(C, p)
    RaptorCodes.r10_hdpc_encode!(C, p)
    return p, d, C
end

function r10_degree()
    p, d, C = init()
    deg = mean(RaptorCodes.trip(x, p)[1] for x in 1:10000)
    println(deg)
end

# RQ degree distribution
function rq_degree()
    p, d, C = init()
    for i in 1:2000
        s = RaptorCodes.lt_generate(C, i, p)
        RaptorCodes.add!(d, s)
    end

    index = Array(1:30)
    probability = [0.005, 0.5, 0.1666, 0.0833, 0.05, 0.0333,
                   0.0238, 0.0179, 0.0139, 0.0111, 0.0091, 0.0076,
                   0.0064, 0.0055, 0.0048, 0.0042, 0.0037, 0.0033,
                   0.0029, 0.0026, 0.0024, 0.0022, 0.002, 0.0018,
                   0.0017, 0.0015, 0.0014, 0.0013, 0.0012, 0.0295]

    println("RQ mean LT degree {}", mean(transpose(index)*probability))

    mean_degree = mean(RaptorCodes.active_degree(row) for row in d.rows)
    println("mean degree is $mean_degree")
end

# K = [4000, 6000, 14000, 34000, 54000, 84000, 134000]
# for k in K
#     p = R10Parameters(k)
#     println("K=$k => $p")
# end

# r10_degree()


# make sure we have a random seed
srand()

function ltqsim(K::Int, overheads::Vector{Int})
    mode = 3998
    delta = 0.9999999701976676
    # delta = 0.999

    # mode = 998
    # delta = 0.9999999701976676
    dd = RaptorCodes.Soliton(K, mode, delta)
    c = LTQ{GF256}(K, dd)
    for _ in 1:100
        df = RaptorCodes.linearsims(overheads, c, 1)
        println(df)
    end
end

function ltqrsim(K::Int, overheads::Vector{Int})
    @assert K == 8000
    # mode = 3998
    mode = K-2
    delta = 0.9999999701976676
    dd = RaptorCodes.Soliton(K, mode, delta)
    c = LTQ{Float64}(K, dd)
    for _ in 1:100
        df = RaptorCodes.linearsims(overheads, c, 100)
        println(df)
    end
end

function ltsim(K::Int, overheads::Vector{Int})
    # mode = 3998
    # delta = 0.9999999701976676
    # delta = 0.999

    mode = 998
    delta = 0.9999999701976676
    dd = RaptorCodes.Soliton(K, mode, delta)
    c = LT(K, dd)
    for _ in 1:100
        df = RaptorCodes.linearsims(overheads, c, 10)
        println(df)
    end
end

function r10sim(K::Int, overheads::Vector{Int})
    c = R10(K)
    for _ in 1:100
        df = RaptorCodes.linearsims(overheads, c, 10)
        println(df)
    end
end

function ldpc10sim(erasures::Vector{Int})
    H = Matrix{Bool}(readdlm("./test/H_612_1224.txt"))
    c = LDPC10{GF256}(H)
    for _ in 1:100
        df = RaptorCodes.linearsims(erasures, c, 1)
        println(df)
    end
end

K = 1224
# reloverheads = linspace(0.3, 0.5, 10)
reloverheads = linspace(0, 0.5, 10)
num_erasures = Vector{Int}(round.(K*reloverheads))
ldpc10sim(num_erasures)

# K = 8000
# reloverheads = linspace(0.25, 0.4, 10)
# overheads = Vector{Int}(round.(K*reloverheads))
# ltqrsim(K, overheads)

# ltsim(K, overheads)
# reloverheads = linspace(0.05, 0.5, 100)
# overheads = Vector{Int}(round.(K*reloverheads))

# overheads = collect(0:20)

# reloverheads = linspace(0, 0.5, 100)
# for K in [20000]
#     overheads = Vector{Int}(round.(K*reloverheads))
#     r10sim(K, overheads)
# end
