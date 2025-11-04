using CairoMakie
using PairPlots

# The simplest table format is just a named tuple of vectors.
# You can also pass a DataFrame, or any other Tables.jl compatible object.
table = (;
    x = randn(10000),
    y = randn(10000),
)

pairplot(table)

##====================================================
using CairoMakie
using PairPlots
# Columns are treated as variables, and rows as samples.
mat = randn(10000,6)
pairplot(mat)

##====================================================

using CairoMakie
using PairPlots

# The simplest table format is just a named tuple of vectors.
# You can also pass a DataFrame, or any other Tables.jl compatible object.
table1 = (;
    x = randn(10000),
    y = randn(10000),
)

table2 = (;
    x = 1 .+ randn(10000),
    y = 2 .+ randn(10000),
    z = randn(10000),
)

fig = pairplot(table1, table2)

##=====================================================

save("myfigure.png", fig)


##=====================================================

using CairoMakie
using PairPlots
using DataFrames

N = 100_000
α = [2randn(N÷2) .+ 6; randn(N÷2)]
β = [3randn(N÷2); 2randn(N÷2)]
γ = randn(N)
δ = β .+ 0.6randn(N)

df = DataFrame(;α, β, γ, δ)

pairplot(df)

pairplot(df, fullgrid=true)

pairplot(df, bottomleft=false, topright=true)

##====================================================

pairplot(df, labels = Dict(
    # basic string
    :α => "parameter 1",
    # Makie rich text
    :β => rich("parameter 2", font=:bold, color=:blue),
    # LaTeX String
    :γ => L"\frac{a}{b}",
))
