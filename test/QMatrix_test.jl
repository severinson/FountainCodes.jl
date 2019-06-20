using FountainCodes, Test

@test (M = QMatrix{Int}(64, 10); M[1,1] == 0)
@test (M = QMatrix{Int}(64, 10); M[1,1] = 0; M[1,1] == 0)
@test (M = QMatrix{Int}(64, 10); M[1,1] = 10; M[1,1] == 10)
@test (M = QMatrix{Int}(64, 10); M[1,1] = 1; M[1,1] == 1)
@test (M = QMatrix{Int}(64, 10); M[1,1] = 1; M[2,2] = 1; subtract!(M,1,2); M[2,1] == 1)
@test (M = QMatrix{Int}(64, 10); M[1,1] = 1; M[2,2] = 1; subtract!(M,5,1,2); M[2,1] == -5)
@test (M = QMatrix{Int}(64, 10); M[1,1] = 1; M[2,1] = 1; FountainCodes.countnz(M,1) == 2)
@test (M = QMatrix{Int}(64, 10); M[64,1] = 1; resize!(M, 128, 10); M[64,1] == 1 && M[65,1] == 0)
@test (M = QMatrix{Int}(64, 10); M[64,10] = 1; resize!(M, 64, 20); M[64,10] == 1 && M[64,11] == 0)
@test (M = QMatrix{Int}(64, 10); M[64,10] = 1; resize!(M, 128, 10); M[64,10] == 1)
@test (M = QMatrix{Int}(64, 10); M[64,10] = 1; resize!(M, 256, 20); M[64,10] == 1)
@test (M = QMatrix{Int}(64, 10); M[64,10] = 2; resize!(M, 256, 20); M[64,10] == 2 && M[65,11] == 0)
@test (M = QMatrix{Int}(64, 10); M[64,10] = 1; M[:,10] = 0; M[64,10] == 0)
@test (M = QMatrix{Int}(64, 10); M[64,10] = 1; M[:,10] = 2; M[1,10] == 2)
@test (M = QMatrix{Int}(64, 10); M[:,1] = 1; M[:,10] = 1; subtract!(M, 1, 1, 10); iszero(getcolumn(M, 1)))

# tests for subtracting one column from another

"""test both cols binary"""
function test_subtract_1()
    m, n = 64, 10
    M = QMatrix{Float64}(m, n)
    M[:, 2] .= true
    subtract!(M, 1, 2)
    if M[:, 1] != trues(m)
        error("expected $(trues(m)), but got $(M[:, 1])")
    end
    return true
end
@test test_subtract_1()

"""test target col binary and source col non-binary"""
function test_subtract_2()
    m, n = 64, 10
    v = 1.2
    M = QMatrix{Float64}(m, n)
    M[:, 2] .= v
    subtract!(M, 1, 2)
    if !all(M[:, 1] .== -v)
        error("expected all values to be $(-v), but got $(M[:, 1])")
    end
    return true
end
@test test_subtract_2()

"""test source col binary and target col non-binary"""
function test_subtract_3()
    m, n = 64, 10
    v = 1.2
    M = QMatrix{Float64}(m, n)
    M[:, 1] .= v
    M[:, 2] .= true
    subtract!(M, 1, 2)
    if !all(M[:, 1] .== (v-1))
        error("expected all values to be $(v-1), but got $(M[:, 1])")
    end
    return true
end
@test test_subtract_3()

"""test both target and source col non-binary"""
function test_subtract_4()
    m, n = 64, 10
    v1 = 1.2
    v2 = 2.0
    M = QMatrix{Float64}(m, n)
    M[:, 1] .= v1
    M[:, 2] .= v2
    subtract!(M, 1, 2)
    if !all(M[:, 1] .== (v1-v2))
        error("expected all values to be $(v1-v2), but got $(M[:, 1])")
    end
    return true
end
@test test_subtract_4()

# tests for subtracting a multiple of one column from another

"""test both cols binary and coef = 1"""
function test_subtract_coef_1()
    m, n = 64, 10
    M = QMatrix{Float64}(m, n)
    M[:, 2] .= true
    subtract!(M, 1.0, 1, 2)
    if !FountainCodes.isbinary(M, 1)
        error("expected column 1 to be binary, but got $(M[:, 1])")
    end
    if M[:, 1] != trues(m)
        error("expected $(trues(m)), but got $(M[:, 1])")
    end
    return true
end
@test test_subtract_coef_1()

"""test both cols binary and coef != 1"""
function test_subtract_coef_2()
    m, n = 64, 10
    M = QMatrix{Float64}(m, n)
    M[:, 2] .= true
    subtract!(M, 2.0, 1, 2)
    if FountainCodes.isbinary(M, 1)
        error("expected column 1 to be non-binary, but got $(M[:, 1])")
    end
    if !all(M[:, 1] .== (-2.0))
        error("expected all values to be -2, but got $(M[:, 1])")
    end
    return true
end
@test test_subtract_coef_2()

"""test target col binary and source col non-binary"""
function test_subtract_coef_3()
    m, n = 64, 10
    v = 1.2
    M = QMatrix{Float64}(m, n)
    M[:, 2] .= v
    subtract!(M, 2.0, 1, 2)
    if !all(M[:, 1] .== -2.4)
        error("expected all values to be -2.4, but got $(M[:, 1])")
    end
    return true
end
@test test_subtract_coef_3()

"""test source col binary and target col non-binary"""
function test_subtract_coef_4()
    m, n = 64, 10
    v = 1.2
    M = QMatrix{Float64}(m, n)
    M[:, 1] .= v
    M[:, 2] .= true
    subtract!(M, 2.0, 1, 2)
    if !all(M[:, 1] .== (1.2-2.0))
        error("expected all values to be $(1.2-2.0), but got $(M[:, 1])")
    end
    return true
end
@test test_subtract_coef_4()

"""test both target and source col non-binary"""
function test_subtract_coef_5()
    m, n = 64, 10
    v1 = 1.2
    v2 = 2.0
    M = QMatrix{Float64}(m, n)
    M[:, 1] .= v1
    M[:, 2] .= v2
    subtract!(M, 10.0, 1, 2)
    if !all(M[:, 1] .== (v1-10*v2))
        error("expected all values to be $(v1-10*v2), but got $(M[:, 1])")
    end
    return true
end
@test test_subtract_coef_5()

# resize tests

"""test that subtracted values persist after resize"""
function test_resize_1()
    M = QMatrix{UInt8}(64, 10)
    M[:,1] = UInt8(1)
    for i in 2:10
        subtract!(M, UInt8(1), i, 1)
    end
    for i in 2:10
        for (j, v) in enumerate(getcolumn(M, i))
            if v != one(v)
                error("M[$j, $i] should be one but is $v")
            end
        end
    end
    resize!(M, 128, 10)
    for i in 2:10
        subtract!(M, UInt8(1), i, 1)
    end
    for i in 2:10
        for (j, v) in enumerate(getcolumn(M, i))
            if !iszero(v)
                error("M[$j, $i] should be zero but is $v")
            end
        end
    end
    return true
end
@test test_resize_1()
