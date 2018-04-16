using RaptorCodes, Base.Test

function test_ldpc_decoder_1()
    H = zeros(Bool, 3, 6)
    H[1, 1] = true
    H[1, 2] = true
    H[2, 3] = true
    H[2, 4] = true
    H[3, 5] = true
    H[3, 6] = true
    println(H)
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


