import RaptorCodes, DataStructures, DataFrames, CSV
@everywhere using RaptorCodes, DataStructures, DataFrames, CSV
@everywhere include("libsim.jl")

""" Packages
- DataFrames
- CSV
"""

function r10_main()
    dir = "results"
    futures = Vector{Future}()
    srs = Vector{Simulation}()
    for K in [1000, 2000, 4000, 8000]
        p = R10Parameters(K)
        C = init(p)
        for overhead in linspace(K*0.05, K*0.5, 10)
            overhead = Int(round(overhead))
            sr = Simulation(p, overhead, C, 100)
            println(repr(sr))
            dct = Dict()
            dct["K"] = K
            dct["overhead"] = overhead
            f = @spawn simulate(sr, dct, dir)
            push!(futures, f)
            push!(srs, sr)
        end
    end
    mkpath(dir)
    results = [fetch(f) for f in futures]
    # for (s, sr) in zip([fetch(f) for f in futures], srs)
    #     filename = joinpath(
    #         dir,
    #         "$(repr(sr)).csv"
    #     )
    #     df = DataFrame(s)
    #     # CSV.write(filename, df)
    # end
end

function lt_main()
    dir = "results"
    futures = Vector{Future}()
    srs = Vector{Simulation}()
    for K in [1000, 2000, 4000, 8000]
        mode = K-1
        delta = 0.0001
        s = RaptorCodes.Soliton(K, mode, delta)
        p = LTParameters(K, s)
        C = init(p)
        for overhead in linspace(K*0.1, K*0.5, 10)
            overhead = Int(round(overhead))
            sr = Simulation(p, overhead, C, 100)
            println(repr(sr))
            dct = Dict()
            dct["K"] = K
            dct["overhead"] = overhead
            dct["mode"] = mode
            dct["delta"] = delta
            f = @spawn simulate(sr, dct, dir)
            push!(futures, f)
            push!(srs, sr)
        end
    end
    mkpath(dir)
    results = [fetch(f) for f in futures]
    # for (s, sr) in zip([fetch(f) for f in futures], srs)
    #     filename = joinpath(
    #         dir,
    #         "$(repr(sr)).csv"
    #     )
    #     df = DataFrame(s)
    #     # CSV.write(filename, df)
    # end
end

r10_main()
