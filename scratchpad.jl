include("simulate.jl")

""" Packages
- DataFrames
- CSV
"""

function r10_main()
    for K in 1000:500:8000
        p = R10Parameters(K)
        C = init(p)
        for overhead in 0:50:500
            sr = Simulation(p, overhead, C, 100)
            println(repr(sr))
            df = simulate(sr, "results")
        end
    end
end

function lt_main()
    futures = Array{Future,1}()
    for K in 1000:500:8000
        mode = K-1
        delta = 0.0001
        s = RaptorCodes.Soliton(K, mode, delta)
        p = LTParameters(K, s)
        C = init(p)
        for overhead in 0:50:500
            sr = Simulation(p, overhead, C, 100)
            println(repr(sr))
            # f = @spawn simulate(sr, "results")
            dct = Dict()
            dct["K"] = K
            dct["overhead"] = overhead
            dct["mode"] = mode
            dct["delta"] = delta
            df = simulate(sr, dct, "results2")
            # push!(futures, f)
        end
    end
    # dfs = Array{DataFrame,1}()
    # for f in futures
    #     df = fetch
    #     push!(dfs, fetch(f))
    # end
end

lt_main()
