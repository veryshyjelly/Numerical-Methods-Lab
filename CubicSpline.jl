function CubicSplineInterploation(x::Vector{Float64}, y::Vector{Float64}, n::Int32)::Vector{Float64}
    # Define the matrix A
    A = zeros(n, n)
    # Define the matrix B
    B = zeros(n)

    # Define functions a, b, c, d to make the code more readable
    a(i) = 4*i 
    b(i) = 4*i + 1
    c(i) = 4*i + 2
    d(i) = 4*i + 3

    # f(xᵢ) = yᵢ
    for i in 1:n-1
        A[a(i), i] = x[i]^3
        A[b(i), i] = x[i]^2
        A[c(i), i] = x[i]
        A[d(i), i] = 1
        B[i] = y[i]
    end
    
    # f(xᵢ₊₁) = yᵢ₊₁
    for i in 1:n-1
        A[a(j), n-1+i] = x[i+1]^3
        A[b(j), n-1+i] = x[i+1]^2
        A[c(j), n-1+i] = x[i+1]
        A[d(j), n-1+i] = 1
        B[n-1+i] = y[i+1]
    end

    # f′(xᵢ) = f′(xᵢ₊₁)
    for i in 1:n-2
        A[a(i), 2*n-2+i] = 3*x[i+1]^2
        A[b(i), 2*n-2+i] = 2*x[i+1]
        A[c(i), 2*n-2+i] = 1
        A[a(i+1), 2*n-2+i] = -3*x[i+1]^2
        A[b(i+1), 2*n-2+i] = -2*x[i+1]
        A[c(i+1), 2*n-2+i] = -1
    end

    # f″(xᵢ) = f″(xᵢ₊₁)
    for i in 1:n-2
        A[a(i), 3*n-4+i] = 6*x[i+1]
        A[b(i), 3*n-4+i] = 2
        A[a(i+1), 3*n-4+i] = -6*x[i+1]
        A[b(i+1), 3*n-4+i] = -2
    end
    
    # f″(x₀) = 0
    A[a(1), 4*n-5] = 6*x[1]
    A[b(1), 4*n-5] = 2
    # f″(xₙ) = 0
    A[a(n-1), 4*n-4] = 6*x[n]
    A[b(n-1), 4*n-4] = 2

    # Solve the system of linear equations
    C = A\B

    # Return the coefficients
    return C

end



