using Plots

# Lex Friedrichs Method
function LFM(u0::Function, f::Function, ΔtByΔx::Float64, x::Vector{Float64}, t::Vector{Float64})::Vector{Float64}
    U::Dict{Float64, Vector{Float64}} = Dict{Float64, Vector{Float64}}()

    merge!(U, Dict(t[1] => u0.(x)))

    # Lax Friedrichs Method
    u(n, j) = (U[t[n-1]][j+1] + U[t[n-1]][j-1])/2 - ΔtByΔx*(f(U[t[n-1]][j+1]) - f(U[t[n-1]][j-1]))/2

    for n in 2:length(t)
        merge!(U, Dict(t[n] => [0; [u(n, j) for j in 2:length(x)-1]; 0]))
        U[t[n]][end] = U[t[n]][end-1]
        U[t[n]][1] = U[t[n]][end]
    end
    
    return U[t[end]]
end

# Define the initial function
function u₀(x::Float64)::Float64
    if abs(x) < (1/3)
        return 1
    else
        return -1
    end
end

# Inviscid Burgers' Equation
f(u) = u^2 / 2

# Define the domain of x
x₀ = -1
xₙ = 1 

# Given value of Δt/Δx 
ΔtByΔx = 0.8

# Define the domain of t
t₀ = 0
tₙ = 0.3

# Solve the PDE for 200 points
Δx = 0.01
Δt = ΔtByΔx*Δx
plot(Vector(x₀:Δx:xₙ), LFM(u₀, f, ΔtByΔx, Vector(x₀:Δx:xₙ), Vector(t₀:Δt:tₙ)),
 label="n = 200", xlims=(x₀, xₙ), ylims=(-1, 1.5), xlabel="x", ylabel="u(x,t)", title="Inviscid Burgers' Equation", color=:blue)


# Solve the PDE for 5000 points
Δx = 0.0004
Δt = ΔtByΔx*Δx
plot!(Vector(x₀:Δx:xₙ), LFM(u₀, f, ΔtByΔx, Vector(x₀:Δx:xₙ), Vector(t₀:Δt:tₙ)),
 label="n = 5000", color=:red)

png("solution2.png")