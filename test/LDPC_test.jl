using RaptorCodes, Base.Test

function test_ldpc_decoder_1()
    H = zeros(Bool, 3, 6)
    H[1, 1] = true
    H[1, 2] = true
    H[2, 3] = true
    H[2, 4] = true
    H[3, 5] = true
    H[3, 6] = true
    c = LDPC10{GF256}(H)
    c.erased[1] = true
    c.erased[3] = true
    c.erased[5] = true
    c.Y[4] = 4
    c.Y[6] = 6
    d = Decoder(c)
    decode!(d)
    get_source!(view(c.Y, c.erased), d)
    if c.Y != Vector{GF256}([0, 0, 4, 4, 6, 6])
        error("LDPC decoding failed")
    end
    return true
end
@test test_ldpc_decoder_1()

function test_ldpc_decoder_2()
    H = Matrix{Bool}(readdlm("./test/H_612_1224.txt"))
    c = LDPC10{GF256}(H)
    c.erased[1:450] = true
    shuffle!(c.erased)
    c.Y[c.erased] = 1 # make sure these get reset to 0 during decoding
    d = Decoder(c)
    decode!(d)
    get_source!(view(c.Y, c.erased), d)
    if !iszero(c.Y)
        error("LDPC decoding failed")
    end
    return true
end
@test test_ldpc_decoder_2()

function test_ldpc_decoder_3()
    H = Matrix{Bool}(readdlm("./test/H_2400_4800.txt"))
    c = LDPC10{GF256}(H)
    c.erased[1:1776] = true
    shuffle!(c.erased)
    c.Y[c.erased] = 1 # make sure these get reset to 0 during decoding
    d = Decoder(c)
    decode!(d)
    get_source!(view(c.Y, c.erased), d)
    if !iszero(c.Y)
        error("LDPC decoding failed")
    end
    return true
end
@test test_ldpc_decoder_3()
