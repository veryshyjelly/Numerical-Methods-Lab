using Plots

function GaussElimination(A::Matrix{Float64}, b::Vector{Float64}, n::Int64)::Vector{Float64}
    for i in 1:n-1
        for j in i+1:n
            m = A[j,i]/A[i,i]
            A[j,:] -= m * A[i,:]
            b[j] -= m * b[i]
        end
    end

    x = zeros(n)
    for i in n:-1:1
        x[i] = (b[i] - dot(A[i,i+1:n], x[i+1:n])) / A[i,i]
    end

    return vec(x)
end

function dot(A::Vector{Float64}, B::Vector{Float64})::Float64
    sum = 0.0
    for i in 1:length(A)
        sum += A[i] * B[i]
    end
    return sum
end

function LeastSquareAppr(x::Vector{Float64}, y::Vector{Float64},n::Int64)::Vector{Float64}
    m = length(x)
    A = zeros(n, n)
    b = zeros(n)

    for i in 1:n
        for j in 1:n
            for k in 1:m
                A[i,j] += x[k]^(2*n-i-j)
            end
        end
        for k in 1:m
            b[i] += y[k] * (x[k]^(n-i))
        end
    end

    return GaussElimination(A, b, n)
end

function Polynomial(coeff::Vector{Float64})::Function
    return x -> sum([coeff[i] * x^(length(coeff)-i) for i in 1:length(coeff)])
end

# Ax::Vector{Float64} = vec([1, 1.5, 2, 2.5, 3, 3.5, 4])
# Ay::Vector{Float64} = vec([25, 31, 27, 28, 36, 35, 32])
# Aresult = LeastSquareAppr(Ax, Ay, 2)
# println("1. y = ax + b, [a, b] = ", Aresult)
# plot(1:0.5:4, Polynomial(Aresult).(1:0.5:4), label="y = ax + b", xlims=(0, 5), ylims=(20, 40))
# scatter!(Ax, Ay, label="data")
# png("1.png")

# Bx::Vector{Float64} = vec([0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
# By::Vector{Float64} = vec([0.0, 0.20, 0.27, 0.30, 0.32, 0.33])
# Bresult = LeastSquareAppr(Bx, By, 4)
# println("2. y = ax³ + bx² + cx + d, [a, b, c, d] = ", Bresult)
# plot(0:0.1:2.5, Polynomial(Bresult).(0:0.1:2.5), label="y = ax³ + bx² + cx + d", xlims=(-0.5, 3), ylims=(0.0, 0.4))
# scatter!(Bx, By, label="data")
# png("2.png")

Sx = Vector(-3:0.4:3)
Sy = sin.(Sx)
SinResult = LeastSquareAppr(Sx, Sy, 10)
# println("3. y = ax⁴ + bx³ + cx² + dx + e, [a, b, c, d, e] = ", SinResult)
plot(-4:0.1:4, Polynomial(SinResult).(-4:0.1:4), label="y = sin(x)", xlims=(-4, 4), ylims=(-1.5, 1.5), color=:red)
SinResult = LeastSquareAppr(Sx, Sy, 4)
plot!(-4:0.1:4, Polynomial(SinResult).(-4:0.1:4), label="y = sin(x)", xlims=(-4, 4), ylims=(-1.5, 1.5), color=:blue)

scatter!(Sx, Sy, label="data")
png("3.png")

