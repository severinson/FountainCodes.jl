using RaptorCodes

function test_gray()
    hw = 3
    for i in RaptorCodes.gray(10, hw)
        if RaptorCodes.hamming_weight(i) != hw
            error("hw($i) != $hw")
        end
    end
    return true
end
@test test_gray()
