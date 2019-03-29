using RaptorCodes
using PyPlot

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

function plot_distribution(c::Union{LT,LTQ})
    trips = [RaptorCodes.trip(X, c) for X in 1:100000]
    samples = [t[1] for t in trips]
    plt[:hist](samples,100)
    # x = linspace(0,2*pi,1000); y = sin.(3*x + 4*cos.(2*x))
    # plot(x, y, color="red", linewidth=2.0, linestyle="--")
    plt[:show]()
end
function plot_symbol_distribution(c::Code)
    samples = Vector{Int}()
    C = [Vector{GF256}([0]) for _ in 1:c.L]
    for X in 1:Int(round(c.L*1.3))
        s = ltgenerate(C, X, c)
        append!(samples, s.neighbours)
    end
    u = unique(samples)
    if length(u) != c.L
        println("$(c.L-length(u)) symbols uncovered")
    end
    # plt[:hist](samples,1000)
    # plt[:show]()
end
function plot_coefficient_distribution(c::LTQ)
    samples = [RaptorCodes.coefficient(X, c) for X in 1:10000]
    plt[:hist](samples,255)
    plt[:show]()
end
# K = 134000
# mode = K-2
# delta = 0.99
# dd = RaptorCodes.Soliton(K, mode, delta)
# c = LTQ{Float64}(K, dd)
# fp = ltfailure_lower(K, 0.3, dd)
# println("fp lower bound = $fp")
# plot_distribution(c)
# plot_symbol_distribution(c)
# plot_coefficient_distribution(c)


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

# function ltqrsim(K::Int, overheads::Vector{Int})
#     mode = K-2
#     delta = 0.999
#     dd = RaptorCodes.Soliton(K, mode, delta)
#     c = LTQ{GF256}(K, dd)
#     for _ in 1:100
#         df = RaptorCodes.linearsims(overheads, c, 1)
#         println(df)
#     end
# end

function ltqrsim(K::Int, mode::Int, delta::Float64, overheads::Vector{Int})
    dd = RaptorCodes.Soliton(K, mode, delta)
    c = LTQ{GF256}(K, dd)
    for _ in 1:100
        df = RaptorCodes.linearsims(overheads, c, 1000)
        println(df)
    end
end

function ltsim(K::Int, overheads::Vector{Int})
    # mode = 3998
    # delta = 0.9999999701976676
    # delta = 0.999

    mode = K-2
    delta = 0.999
    # delta = 0.9999999701976676
    dd = RaptorCodes.Soliton(K, mode, delta)
    c = LT(K, dd)
    for _ in 1:1000
        df = RaptorCodes.linearsims(overheads, c, 100)
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

function RQsim(K::Int, overheads::Vector{Int})
    c = RQ(K)
    for _ in 1:1000
        df = RaptorCodes.linearsims(overheads, c, 100)
        println(df)
    end
end

function ldpc10sim(erasures::Vector{Int})
    H = Matrix{Bool}(readdlm("./test/H_612_1224.txt"))
    # H = Matrix{Bool}(readdlm("./test/H_2400_4800.txt"))
    c = LDPC10{GF256}(H)
    for _ in 1:1000
        df = RaptorCodes.linearsims(erasures, c, 10000)
        println(df)
    end
end

# K = 4800
# reloverheads = linspace(0.3, 0.5, 10)
# num_erasures = Vector{Int}(round.(K*reloverheads))
# num_erasures = [1973, 2080]#, 2187, 2240] #, 2293, 2347, 2400]
# num_erasures = [490, 503]#, 517]#, 530, 544]#, 558, 571]#, 585, 598]
# ldpc10sim(num_erasures)

# K = 8000
# reloverheads = linspace(0.25, 0.4, 10)
# overheads = Vector{Int}(round.(K*reloverheads))
# ltqrsim(K, overheads)

# ltsim(K, overheads)
# reloverheads = linspace(0.05, 0.5, 100)
# overheads = Vector{Int}(round.(K*reloverheads))

# overheads = collect(0:20)

# K = 4015
# reloverheads = linspace(0, 0.1, 10)
# overheads = Vector{Int}(round.(K*reloverheads))
# RQsim(K, overheads)
# println(RQ(1000).P)
# println(RQ(10000).P)
# println(RQ(50000).P)

# for K in [20000]
#     overheads = Vector{Int}(round.(K*reloverheads))
#     r10sim(K, overheads)
# end

## partition plot
# rate 2/3, 30% overhead
# 24, 3000

# K = 6000
# reloverheads = [0.3]
# overheads = Vector{Int}(round.(K*reloverheads))
# ltqrsim(K, 5998, 0.9999999701976676, overheads)

# K = 24
# reloverheads = [0.3]
# overheads = Vector{Int}(round.(K*reloverheads))
# ltqrsim(K, 22, 0.9999999701976676, overheads)

## size plot
# rate 2/3, 30% overhead
# 4000, 6000, 14000, 34000, 54000, 84000, 134000

# num_inputs = [10, 24, 140, 850, 2160, 5250, 13400, 4000, 6000, 14000, 34000, 54000, 84000, 134000]
# num_inputs = [10, 24, 140, 850, 2160]
# num_inputs = [5250, 13400, 4000]
# num_inputs = [6000, 14000, 34000]
# num_inputs = [54000]
# num_inputs = [84000]
# num_inputs = [134000]

# num_inputs = [
#     12600
#     14424
#     15884
#     17080
#     18156
#     21460
#     22176
#     23712
#     25568
#     26600
#     28968
#     32040
#     36750
#     44408
#     59800
#     10
#     24
#     44
#     70
#     102
#     290
#     352
#     494
#     752
#     1064
#     1704
#     2670
#     5250
#     11102
#     29900
# ]

num_inputs = [26, 156, 572, 1430, 2574]
# num_inputs = [50000, 49998, 50076, 49764, 50050, 48906]

# num_inputs = [134000]
# for K in num_inputs
#     reloverheads = [0.335]
#     overheads = Vector{Int}(round.(K*reloverheads))
#     # ltqrsim(K, K-2, 0.9999999701976676, overheads)
#     # ltqrsim(K, K-2, 0.3001587688922852, overheads) # 13400, 1e-3
#     # ltqrsim(K, K-2, 0.533183008432383, overheads) # 134000, 1e-3
#     # ltqrsim(K, K-2, 0.0016897618770599196, overheads) # 13400, 1e-6
#     # ltqrsim(K, K-2, 0.0030159652233123476, overheads) # 134000, 1e-6
#     ltqrsim(K, K-2, 1.7076730728149245e-5, overheads) # 134000, 1e-9
# end

# K = 2
# reloverheads = [0.335]
# overheads = Vector{Int}(round.(K*reloverheads))
# ltqrsim(2, 2, 0.9999999701976676, overheads)

# LT code parameter plot
# K = 2400

# # 1.2/1e-1
# reloverheads = [0.2]
# overheads = Vector{Int}(round.(K*reloverheads))
# delta = 0.9999999701976676
# mode = 2398
# ltqrsim(K, mode, delta, overheads)

# # 1.3/1e-1
# reloverheads = [0.3]
# overheads = Vector{Int}(round.(K*reloverheads))
# delta = 0.9999999701976676
# mode = 2398
# ltqrsim(K, mode, delta, overheads)

# # 1.2/1e-3
# reloverheads = [0.2]
# overheads = Vector{Int}(round.(K*reloverheads))
# delta=0.058706909418105496
# mode = 2398
# ltqrsim(K, mode, delta, overheads)

# # 1.3/1e-3
# reloverheads = [0.3]
# overheads = Vector{Int}(round.(K*reloverheads))
# delta = 0.14514568448066567
# mode = 2398
# ltqrsim(K, mode, delta, overheads)

# K = 2400
# # 1.2/1e-1
# reloverheads = [0.3]
# overheads = Vector{Int}(round.(K*reloverheads))
# delta = 0.9999999701976676
# mode = K-2
# ltqrsim(K, mode, delta, overheads)

# # 1.3/1e-1
# reloverheads = [0.37]
# overheads = Vector{Int}(round.(K*reloverheads))
# delta = 0.9999999701976676
# mode = K-2
# ltqrsim(K, mode, delta, overheads)

# # 1.2/1e-3
# reloverheads = [0.3]
# overheads = Vector{Int}(round.(K*reloverheads))
# delta = 0.14514568448066567
# mode = K-2
# ltqrsim(K, mode, delta, overheads)

# # 1.3/1e-3
# reloverheads = [0.37]
# overheads = Vector{Int}(round.(K*reloverheads))
# delta = 0.2555731832981084
# mode = K-2
# ltqrsim(K, mode, delta, overheads)
