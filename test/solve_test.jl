using RaptorCodes, Nemo

function test_solve_dense_1()
    R, x = Nemo.FiniteField(2, 8, "x")
    num_inputs = 5
    A = zero(MatrixSpace(R, 4, 2))
    v = zero(MatrixSpace(R, 2, 1))
    A[1,1] = x^3+x^2
    A[2,1] = x^3+x^2
    A[2,2] = x^6+x^5+x^3+1
    A[3,2] = x^6+x^5+x^4+x^3+x^2+x
    A[4,1] = x^2
    A[4,2] = x^5
    v[1, 1] = x
    v[2, 1] = x^2
    b = A*v
    vv = RaptorCodes.solve(A, b)
    if A*vv != b
        error("incorrect solution to system")
    end
    return true
end
@test test_solve_dense_1()

function test_solve_dense_2()
    R, x = Nemo.FiniteField(2, 8, "x")
    num_inputs = 5
    A = zero(MatrixSpace(R, 4, 2))
    v = zero(MatrixSpace(R, 2, 1))
    A[1, 2] = 1
    A[4, 1] = 1
    v[1, 1] = x
    v[2, 1] = x^2
    b = A*v
    vv = RaptorCodes.solve(A, b)
    if A*vv != b
        error("incorrect solution to system")
    end
    return true
end
@test test_solve_dense_2()

function test_solve_dense_3()
    R, x = Nemo.FiniteField(2, 8, "x")
    num_inputs = 5
    A = zero(MatrixSpace(R, 2, 2))
    v = zero(MatrixSpace(R, 2, 1))
    A[1, 1] = 1
    A[1, 2] = x
    A[2, 1] = x+x^3
    A[2, 2] = x^4
    v[1, 1] = x
    v[2, 1] = x^2
    b = A*v
    vv = RaptorCodes.solve(A, b)
    if A*vv != b
        error("incorrect solution to system")
    end
    return true
end
@test test_solve_dense_3()
