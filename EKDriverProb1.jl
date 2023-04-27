using SpecialFunctions, Printf, JLD2, CSV, DataFrames

include("functionsProb1.jl")

# Specify parameter values
# ------------------------------------------------------------------
N = 3

L = ones(N)
A = Array{Float64, 1}(undef, N)
d = Array{Float64, 2}(undef, N,N)

for n = 1:N
    A[n] = n
    for i in 1:N
        if n == i
            d[n,i] = 1.0
        else
            d[n,i] = 1.0
        end # end if
    end # end for i in 1:N
end # end for n in 1:N

theta = 4.0
eta = 2.0
v = 0.4

gam = gamma(1.0+(1.0.-eta)./theta).^(1.0./(1.0.-eta))
# ------------------------------------------------------------------



# Iterate on wages
# ------------------------------------------------------------------

# Initial guess for wage vector
w0 = rand(N)
w0 = w0./sum(w0.*L) # Normalize s.t. world GDP=1
P0 = rand(N)

# Compute equilibrium wage
@time P, bts, u, w, EX, IM = compute_wage(w0::Array{Float64,1},
                                A::Array{Float64,1},
                                d::Array{Float64,2},
                                theta::Float64,
                                gam::Float64,
                                N::Int64,
                                v::Float64,
                                P0::Array{Float64,1})

print(w./P)
# ------------------------------------------------------------------




#@save("model_parms.jld2", N, L, A, d, theta, eta, gam, v)

# @save("model_gen_data.jld2", w, P, bts)
#
#
# CSV.write("bilateral_trade_shares.csv", DataFrame(reshape(bts,N*N,1)),
#                       header = false)
#
# CSV.write("prices.csv", DataFrame(reshape(P,N,1)),
#                          header = false)
#
# CSV.write("wages.csv", DataFrame(reshape(w,N,1)),
#                          header = false)
