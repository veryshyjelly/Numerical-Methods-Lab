function CubicSplineInterploation(x::Vector{Float64}, y::Vector{Float64}, n::Int64)::Vector{Float64}
    # Define the matrix A
    A = zeros(4*n-4, 4*n-4)
    # Define the matrix B
    B = zeros(4*n-4)

    # Define functions a, b, c, d to make the code more readable
    a(i) = 4*(i-1) + 1 
    b(i) = 4*(i-1) + 2
    c(i) = 4*(i-1) + 3
    d(i) = 4*(i-1) + 4

    # f(xᵢ) = yᵢ
    for i in 1:n-1
        A[i, a(i)] = x[i]^3
        A[i, b(i)] = x[i]^2
        A[i, c(i)] = x[i]
        A[i, d(i)] = 1
        B[i] = y[i]
    end
    
    # f(xᵢ₊₁) = yᵢ₊₁
    for i in 1:n-1
        A[n-1+i, a(i)] = x[i+1]^3
        A[n-1+i, b(i)] = x[i+1]^2
        A[n-1+i, c(i)] = x[i+1]
        A[n-1+i, d(i)] = 1
        B[n-1+i] = y[i+1]
    end

    # f′(xᵢ) = f′(xᵢ₊₁)
    for i in 1:n-2
        A[2*n-2+i, a(i)] = 3*x[i+1]^2
        A[2*n-2+i, b(i)] = 2*x[i+1]
        A[2*n-2+i, c(i)] = 1
        A[2*n-2+i, a(i+1)] = -3*x[i+1]^2
        A[2*n-2+i, b(i+1)] = -2*x[i+1]
        A[2*n-2+i, c(i+1)] = -1
    end

    # f″(xᵢ) = f″(xᵢ₊₁)
    for i in 1:n-2
        A[3*n-3+i, a(i)] = 6*x[i+1]
        A[3*n-3+i, b(i)] = 2
        A[3*n-3+i, a(i+1)] = -6*x[i+1]
        A[3*n-3+i, b(i+1)] = -2
    end
    
    # f″(x₀) = 0
    A[4*n-5, 1] = 6*x[1]
    A[4*n-5, b(1)] = 2
    # f″(xₙ) = 0
    A[4*n-4, a(n-1)] = 6*x[n]
    A[4*n-4, b(n-1)] = 2

    # Solve the system of linear equations
    show(stdout, "text/plain", A)
    show(stdout, "text/plain", B)
    C = A\B

    # Return the coefficients
    return C

end

CubicPol(a, b, c, d) = x -> a*x^3 + b*x^2 + c*x + d

f(x) = x^3 + 2
X::Vector{Float64} = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
Y::Vector{Float64} = f.(X)

C = CubicSplineInterploation(X, Y, 6)

f₁ = CubicPol(C[1], C[2], C[3], C[4])
f₂ = CubicPol(C[5], C[6], C[7], C[8])
f₃ = CubicPol(C[9], C[10], C[11], C[12])
f₄ = CubicPol(C[13], C[14], C[15], C[16])
f₅ = CubicPol(C[17], C[18], C[19], C[20])

print(f(0.1) - f₁(0.1))
print(f(0.3) - f₂(0.3))
print(f(0.5) - f₃(0.5))

println(C)
