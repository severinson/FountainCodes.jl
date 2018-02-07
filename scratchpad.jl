using RaptorCodes, DataStructures, DataFrames, CSV

""" Packages
- DataFrames
- CSV
"""

struct Simulation
    p::RaptorCodes.Parameters
    overhead::Int # absolute reception overhead
    C::Array{RaptorCodes.ISymbol,1}
    n::Int # number of samples
end

function repr(r::Simulation)
    return "Simulation($(RaptorCodes.repr(r.p)), $(r.overhead))"
end

function init(p::R10Parameters)
    C = Array{RaptorCodes.ISymbol,1}(p.L)
    for i = 1:p.K
        C[i] = RaptorCodes.ISymbol(i, Set([i]))
    end
    RaptorCodes.r10_ldpc_encode!(C, p)
    RaptorCodes.r10_hdpc_encode!(C, p)
    return C
end

function init(p::LTParameters)
    C = Array{RaptorCodes.ISymbol,1}(p.L)
    for i = 1:p.K
        C[i] = RaptorCodes.ISymbol(i, Set([i]))
    end
    return C
end

function sample(sr::Simulation)
    d = RaptorCodes.Decoder(sr.p)
    for _ in 1:(sr.p.K+sr.overhead)
        X = rand(1:(2<<16))
        s = RaptorCodes.lt_generate(sr.C, X, sr.p)
        RaptorCodes.add!(d, s)
    end
    RaptorCodes.decode!(d, false)
    return d.metrics
end

function simulate(sr::Simulation, dir::String)
    mkpath(dir)
    filename = joinpath(
        dir,
        "$(repr(sr)).csv"
    )
    try
        return CSV.read(filename)
    catch
    end
    samples = DefaultDict(Array{Float64,1})
    for i in 1:sr.n
        for v in sample(sr)
            push!(samples[v[1]], v[2])
        end
    end
    df = DataFrame(samples)
    CSV.write(filename, df)
    return df
end

function r10_main()
    for K in 1000:500:8000
        p = R10Parameters(K)
        C = init(p)
        for overhead in 0:50:500
            sr = Simulation(p, overhead, C, 100)
            println(repr(sr))
            df = @spawn simulate(sr, "results")
        end
    end
end

function lt_main()
    for K in 1000:500:8000
        s = RaptorCodes.Soliton(K, K-1, 0.0001)
        p = LTParameters(K, s)
        C = init(p)
        for overhead in 0:50:500
            sr = Simulation(p, overhead, C, 100)
            println(repr(sr))
            df = @spawn simulate(sr, "results")
        end
    end
end

lt_main()
