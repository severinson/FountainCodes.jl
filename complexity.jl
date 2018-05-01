# Compute the encoding complexity of R10/RQ codes.
# Store as a dataframe.

using RaptorCodes, JuliaDB, IterableTables, FileIO

"""
    complexity(c::R10)

Return a tuple (b, a, m) with the number of bit XORs, field additions, and field
multiplications required for the precode.

R10 degree distribution mean: 4.631353378295898

"""
function complexity(c::R10)
    a = 0.0

    # source symbols are connected to 3 LDPC symbols
    a += c.K * 3

    # each HDPC symbol has about 50% non-zero entries
    w = (c.K+c.S)/2
    a += w*c.H

    return a, 0.0, 0.0
end

"""
    complexity(c::RQ)

Return a tuple (b, a, m) with the number of bit XORs, field additions, and field
multiplications required for the precode.

RQ degree distribution mean w/o PI: 4.816641807556152
RQ degree distribution mean: 7.152566

"""
function complexity(c::RQ)
    b, a, m = 0.0, 0.0, 0.0

    # source symbols are connected to 3 LDPC symbols
    b += 3c.K

    # each LDPC symbol is connected to 2 PI symbols
    b += 2c.S

    # each HDPC symbol is composed of H+S GF256 entries
    a += (c.K+c.S)*c.H
    m = a

    return b, a, m
end

"""
    RQavg_empiric()

Compute the average degree empirically. Takes into account the degree
distribution over the PI symbols.

"""
function RQavg_empiric()
    r = 0
    n = 1000000
    c = RQ(1000)
    for X in 1:n
        d, _, _, d1, _, _ = RaptorCodes.RQ_tuple(X, c)
        r += d+d1
    end
    println(r/n)
end

function R10avg()
    dt = [0, 1, 2, 3, 4, 10, 11, 40]
    ft = [0, 10241, 491582, 712794, 831695, 948446, 1032189, 1048576]
    r = 0.0
    for i in 2:length(dt)
        d = dt[i]
        f = ft[i]-ft[i-1]
        r += d*f/ft[end]
    end
    return r
end

function RQavg()
    dt = Vector(0:30)
    ft = [0, 5243, 529531, 704294, 791675, 844104, 879057, 904023, 922747, 937311,
          948962, 958494, 966438, 973160, 978921, 983914, 988283, 992138, 995565,
          998631, 1001391, 1003887, 1006157, 1008229, 1010129, 1011876, 1013490,
          1014983, 1016370, 1017662, 1048576]
    r = 0.0
    for i in 2:length(dt)
        d = dt[i]
        f = ft[i]-ft[i-1]
        r += d*f/ft[end]
    end
    return r
end

function main()
    b = Vector{Float64}(289)
    a = Vector{Float64}(289)
    m = Vector{Float64}(289)
    K = Vector{Float64}(289)
    for (i, Kp) in enumerate(RaptorCodes.RQ_parameter_table[1:289,1])
        c10 = R10(Kp)
        b[i], a[i], m[i] = complexity(c10)
        K[i] = Kp
    end
    t = table(K, b, m, names=[:K, :b, :f])
    FileIO.save("R10.csv", t)
    println(t)

    n = length(RaptorCodes.RQ_parameter_table[:,1])
    b = Vector{Float64}(n)
    a = Vector{Float64}(n)
    m = Vector{Float64}(n)
    K = Vector{Float64}(n)
    for (i, Kp) in enumerate(RaptorCodes.RQ_parameter_table[:,1])
        cQ = RQ(Kp)
        b[i], a[i], m[i] = complexity(cQ)
        K[i] = Kp
    end
    t = table(K, b, m, names=[:K, :b, :f])
    FileIO.save("RQ.csv", t)
    println(t)
end

# RQavg_empiric()
main()

