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

function qltsim(overheads::Vector{Int})
    mode = 3998
    delta = 0.9999999701976676
    dd = RaptorCodes.Soliton(K, mode, delta)
    c = QLTParameters(K, dd)
    for _ in 1:100
        df = RaptorCodes.linearsims(overheads, c, 100)
        println(df)
    end
end

function ltsim(overheads::Vector{Int})
    mode = 3998
    delta = 0.9999999701976676
    dd = RaptorCodes.Soliton(K, mode, delta)
    c = LTParameters(K, dd)
    for _ in 1:100
        df = RaptorCodes.linearsims(overheads, c, 100)
        println(df)
    end
end

function r10sim(overheads::Vector{Int})
    c = R10Parameters(K)
    for _ in 1:100
        df = RaptorCodes.linearsims(overheads, c, 100)
        println(df)
    end
end

K = 4000
a = 0.01
c = 1.5
reloverheads = [a*c^i for i in 0:10]
overheads = Vector{Int}(round.(K*reloverheads))


reloverheads = linspace(0.25, 0.4, 100)
overheads = Vector{Int}(round.(K*reloverheads))
# qltsim(overheads)
# ltsim(overheads)
reloverheads = linspace(0, 0.05, 100)
overheads = Vector{Int}(round.(K*reloverheads))
r10sim(overheads)
