using Nemo

doc"Solve Ax=b for x."
function solve(A, b)
    if rows(A) != rows(b)
        error("A and b must have equal number of rows")
    end
    if cols(b) != 1
        error("b must have exactly 1 column")
    end

    m = rows(A)
    n = cols(A)

    # concatenate the columns of A and b
    M = similar(A, m, n+1)
    for i = 1:m
        for j = 1:n
            M[i, j] = A[i, j]
        end
        M[i, n+1] = b[i, 1]
    end

    # solve the system of equations
    r = rref!(M)
    if r < cols(A)
        error("system of equations is under-determined")
    end
    if r > cols(A)
        error("system of equations is over-determined")
    end

    # copy the solution into x
    x = similar(A, n, 1)
    for i = 1:n
        x[i, 1] = M[i, n+1]
    end
    return x
end
