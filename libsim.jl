struct Simulation{VT<:RaptorCodes.Value}
    p::RaptorCodes.Parameters
    overhead::Int # absolute reception overhead
    C::Vector{RaptorCodes.ISymbol{VT}}
    n::Int # number of samples
end

function repr(r::Simulation)
    return "Simulation($(RaptorCodes.repr(r.p)), $(r.overhead))"
end

function init(p::R10Parameters)
    C = Vector{RaptorCodes.ISymbol{R10Value}}(p.L)
    for i = 1:p.K
        C[i] = RaptorCodes.ISymbol(R10Value(i), Set([i]))
    end
    RaptorCodes.r10_ldpc_encode!(C, p)
    RaptorCodes.r10_hdpc_encode!(C, p)
    return C
end

function init(p::LTParameters)
    C = Vector{RaptorCodes.ISymbol{R10Value}}(p.L)
    for i = 1:p.K
        C[i] = RaptorCodes.ISymbol(R10Value(i), Set([i]))
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

function simulate(sr::Simulation, dct::Dict, dir::String)
    mkpath(dir)
    filename = joinpath(
        dir,
        "$(repr(sr)).csv"
    )
    try
        return CSV.read(filename)
    catch
    end
    println("starting simulation for $(repr(sr))")
    samples = DefaultDict(Vector)
    for i in 1:sr.n
        for v in dct
            push!(samples[v[1]], v[2])
        end
        for v in sample(sr)
            push!(samples[v[1]], v[2])
        end
    end
    df = DataFrame(samples)
    CSV.write(filename, df)
    println("finished simulation for $(repr(sr))")
    return samples
end
