using BenchmarkTools
using DelimitedFiles
using ICS
using LinearAlgebra
using Random

function tcov_original(X::Union{Matrix{Float64}, DataFrame}, beta=2)
 
    if X isa DataFrame
        X = Matrix(X)
    end

    n, p = size(X)
    b = -beta / 2.0

    cov_inv = inv(ICS.cov(X, dims=1))

    V = zeros(p, p)
    denominator = 0.0

    for i in 2:n
        for j in 1:(i-1)
            diff = X[i, :] .- X[j, :]
            r_sq = dot(diff, cov_inv * diff)
            w = exp(b * r_sq)
            V .+= w * (diff * diff')
            denominator += w
        end
    end

    return Scatter(nothing, V/denominator, "TCOV")
end

X = rand(1000, 10)

bb1 = @benchmark tcov(X) samples=19
bb2 = @benchmark tcov_original(X) samples=19

writedlm("tcov-julia.csv", bb1.times, ',')
writedlm("tcov_original-julia.csv", bb2.times, ',')
