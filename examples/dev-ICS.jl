using Revise, Pkg
Pkg.activate(".julia/dev/ICS")
using ICS

## Workflow for parallel development with Robustbase
using Revise, Pkg
Pkg.develop(path=".julia/dev/Robustbase")        # added Robustbase to the development
Pkg.activate(".julia/dev/ICS")                   # added ICS 
using ICS

## To run the tests:
##  Change to the root directory of the package, activate the package and run 'test'
cd("C:/Users/valen/.julia/dev/ICS")

import Pkg; Pkg.activate("."); Pkg.test()

## To reproduce the CI
import Pkg; Pkg.test(;coverage=true, julia_args=["--check-bounds=yes", "--compiled-modules=yes", "--depwarn=yes"], force_latest_compatible_version=false, allow_reresolve=true)

## To build documentation ...
##  Change to the root directory of the package, activate the environment 'docs' and include make.jl
cd("C:/Users/valen/.julia/dev/ICS")
using Pkg
Pkg.activate("docs")
include("C:/Users/valen/.julia/dev/ICS/docs/make.jl")


## To view the built static documentation:
using LiveServer
servedocs()

##============================================================
using RDatasets
iris = dataset("datasets", "iris")
x = iris[:,1:4]
X = Matrix(iris[:,1:4])

## Scatter
X = Robustbase.hbk[:,1:3]

using ICS
cov2(X)
covW(X)
covAxis(X)
cov4(X)
tcov(X)

mcd_raw(X)
mcd_rwt(X)

##  ICSModel
ics = ICSModel(S2=covW)
fit!(ics, x)
scores = predict(ics, X)

ics = ICSModel()
scores = fit_predict!(ics, X)

##=============================================
using DataFrames
using RCall
using Test

## Load R libraries 
R"library(ICSOutlier)"
R"library(ICSClust)"

function doTestScatter(X::Union{Matrix{Float64}, DataFrame}; which::String="cov2")
    if which == "cov2"
        cc=R"ICS_cov($X, location=TRUE)"
        ss=cov2(X)
    elseif which == "cov4"
        cc=R"ICS_cov4($X, location='mean')"
        ss=cov4(X)
    elseif which == "covW"
        cc=R"ICS_covW($X, location=TRUE, alpha=1, cf=1)"
        ss=covW(X, alpha=1, cf=1)
    elseif which == "covAxis"
        cc=R"ICS_covAxis($X, location=TRUE)"
        ss=covAxis(X)
    elseif which == "tcov"
        cc=R"ICS_tcov($X, beta=2)"
        ss=tcov(X, 2.0)
    elseif which == "mcd_raw"
        cc=R"ICS_mcd_raw($X, location=TRUE, nsamp=500, alpha=0.5)"
        ss=mcd_raw(X, nsamp=500, alpha=0.5)
    else
        error("Undefined scatter:", which)       
    end
    
    @test(ss.label == rcopy(cc)[Symbol("label")])
    if which != "tcov" 
        @test(isapprox(ss.location, rcopy(cc)[Symbol("location")]))
    end
    @test(isapprox(ss.scatter, rcopy(cc)[Symbol("scatter")]))
end

X = rcopy(R"data('hbk', package='robustbase'); x=hbk[,1:3]");
doTestScatter(X, which="cov2")
doTestScatter(X, which="cov4")
doTestScatter(X, which="covW")
doTestScatter(X, which="covAxis")
doTestScatter(X, which="tcov")
##  doTestScatter(X, which="mcd_raw")

using StatsPlots
X = rcopy(R"data('hbk', package='robustbase'); x=hbk[,1:3]");
ics = ICSModel();
scores = fit_predict!(ics, X);
corrplot(scores, label=["x$i" for i=1:size(X, 2)])

using CairoMakie
using PairPlots

ics = ICSModel(S2=cov4, algorithm="whiten");
scores = fit_predict!(ics, X);

pairplot(scores)
pairplot((scores[1:14,:], scores[15:75,:]))

pairplot(scores[1:14,:], scores[15:75,:] => 
    (PairPlots.Scatter(color=:orange),
     PairPlots.MarginDensity(),))


ics = ICSModel(algorithm="standard");
scores = fit_predict!(ics, X);


##=============================================
##
##  HTP3 data
##
using DataFrames
using RCall
using Test
using ICS

## Load R libraries 
R"library(ICSOutlier)"
R"library(ICSClust)"

X = rcopy(R"data('HTP3', package='ICSOutlier');x=HTP3");
c2 = cov2(X);
c4 = cov4(X);       # this does not work in R
cw = covW(X);       # this does not work in R

rics_standard=R"try(ics_res_whiten <- ICS($X, algorithm = 'standard'))"
ics_standard = ICSModel(algorithm="standard");
scores = fit_predict!(ics_standard, X);

rics_whiten=R"try(ics_res_whiten <- ICS($X, algorithm = 'whiten'))"
ics_whiten = ICSModel(algorithm="whiten");
scores = fit_predict!(ics_whiten, X);

ric_qr=R"ICS($X, algorithm = 'QR')"
ics_qr = ICSModel(S2=cov4, algorithm="QR");
scores = fit_predict!(ics_qr, X);

@test(isapprox(ics_qr.kurtosis_, rcopy(ric_qr)[Symbol("gen_kurtosis")]))
@test(isapprox(ics_qr.skewness_, rcopy(ric_qr)[Symbol("gen_skewness")]))
@test(isapprox(ics_qr.W_, rcopy(ric_qr)[Symbol("W")]))
@test(isapprox(scores, rcopy(ric_qr)[Symbol("scores")]))


## 1. Scree plot
fig = scree_plot(ics_qr)

using StatsPlots
savefig(fig, "C:/projects/statproj/julia/HTP3_kurtosis.png")


## 2. Outlier detection plot (index plot)
using StatsPlots
clusters = Vector{String}(undef, size(scores,1))
fill!(clusters, "normal")
clusters[32] = "outlier"

fig=outlier_plot(ics_qr; clusters=clusters, legend=:topright)


savefig(fig, "C:/projects/statproj/julia/HTP3_outlier.png")

ric_qr=R"ICS($X, algorithm = 'QR', fix_signs='W')"
ics_qr = ICSModel(S2=cov4, algorithm="QR", fix_signs="W");
scores = fit_predict!(ics_qr, X);

@test(isapprox(ics_qr.kurtosis_, rcopy(ric_qr)[Symbol("gen_kurtosis")]))
## both are Nothing
@test(ics_qr.skewness_ == rcopy(ric_qr)[Symbol("gen_skewness")])
@test(isapprox(ics_qr.W_, rcopy(ric_qr)[Symbol("W")]))
@test(isapprox(scores, rcopy(ric_qr)[Symbol("scores")]))



rii=R"rx=ICS($X, algorithm = 'QR')"
row_signs=R"rx=ICS($X, algorithm = 'QR'); apply(rx$W, 1L, ICS:::.sign.max)"
row_norms=R"rx=ICS($X, algorithm = 'QR'); apply(rx$W, 1L, ICS:::.sign.max); sqrt(rowSums(rx$W^2))"
W_final=R"rx=ICS($X, algorithm = 'QR'); row_signs=apply(rx$W, 1L, ICS:::.sign.max); row_norms=sqrt(rowSums(rx$W^2)); W_final=sweep(rx$W, 1L, row_norms * row_signs, '/')"

ii = ICSModel(S2=cov4, algorithm="QR");
scores = fit_predict!(ii, X);
ii_row_signs=[ICS.sign_max(ii.W_[i,:]) for i in 1:size(ii.W_,1)]
ii_row_norms=sqrt.(sum(ii.W_.^2, dims=2))
ii_W_final = (ii.W_ ./ (ii_row_signs .* ii_row_norms))

@test(isapprox(ii.W_, rcopy(rii)[Symbol("W")]))
@test(isapprox(rcopy(row_signs), ii_row_signs))
@test(isapprox(rcopy(row_norms), ii_row_norms))
@test(isapprox(rcopy(W_final), ii_W_final))

clusters = Vector{String}(undef, size(scores,1))
fill!(clusters, "normal")
clusters[32] = "outlier"

fig=outlier_plot(ics_qr; clusters=clusters, legend=:topright)

##==================================================================
##
##  data_philips from cellWise

using DataFrames
using RCall
using Test
using ICS

## Load R libraries 
R"library(ICSOutlier)"
R"library(ICSClust)"

R"library(cellWise)"
X = rcopy(R"data('data_philips', package='cellWise'); x=data_philips");

rics=R"ICS($X, S1=ICS_mcd_raw, S2=ICS_cov, S1_args = list(alpha = 0.50, nsamp = 500))"
ics = ICSModel(S1=mcd_raw, S2=cov2);
scores = fit_predict!(ics, X);

#=
@test(isapprox(ics.kurtosis_, rcopy(rics)[Symbol("gen_kurtosis")]))
@test(isapprox(ics.skewness_, rcopy(rics)[Symbol("gen_skewness")]))
@test(isapprox(ics.W_, rcopy(rics)[Symbol("W")]))
@test(isapprox(scores, rcopy(rics)[Symbol("scores")]))
=#

scree_plot(ics)

clusters = repeat(["Group1"], outer=size(scores,1))
clusters[1:100] .= "Group2"
clusters[491:565] .= "Group3"

component_plot2(ics, clusters=clusters)
component_plot2(ics, clusters=clusters, select=[4])


## data_philips with center=TRUE

rics=R"ICS($X, S1=ICS_mcd_raw, S2=ICS_cov, S1_args = list(alpha = 0.50, nsamp = 500), 
    center = TRUE, algorithm = 'standard')";
ics = ICSModel(S1=mcd_raw, S2=cov2, center=true, algorithm="standard", 
    S1_args=Dict{Symbol, Any}(:alpha=>0.5, :nsamp=>500));
scores = fit_predict!(ics, X);

#=
@test(isapprox(ics.kurtosis_, rcopy(rics)[Symbol("gen_kurtosis")]))
@test(isapprox(ics.skewness_, rcopy(rics)[Symbol("gen_skewness")]))
@test(isapprox(ics.W_, rcopy(rics)[Symbol("W")]))
@test(isapprox(scores, rcopy(rics)[Symbol("scores")]))
=#

scree_plot(ics)
component_plot2(ics, clusters=clusters)



using CairoMakie
using PairPlots
select = 1:2
ss = DataFrame(scores, ["IC$i" for i=1:size(scores,2)])

pairplot(scores)
fig=pairplot(ss[clusters .== "Group1", select], 
          ss[clusters .== "Group2", select],
          ss[clusters .== "Group3", select] =>
         (PairPlots.Scatter(),
          PairPlots.MarginDensity(),), 
          fullgrid=true)

pairplot(scores => (PairPlots.Scatter(markersize=8, color=:red, group=clusters),
                    PairPlots.MarginDensity(),))

save("Philips.png", fig)

select = vcat(1:3, 7:9)
fig=pairplot(ss[clusters .== "Group1", select], 
          ss[clusters .== "Group2", select],
          ss[clusters .== "Group3", select] =>
         (PairPlots.Scatter(),
          PairPlots.MarginDensity(),), 
          fullgrid=true)

save("Philips-all.png", fig)

##==================================================================
##
##  penguins

using DataFrames
using RCall
using Test

## Load R libraries 
R"library(ICSOutlier)"
R"library(ICSClust)"

R"library(palmerpenguins)"
penguins = X = rcopy(R"data('penguins', package='palmerpenguins'); x=penguins; x=x[-c(4,9,10,11,12,48,179,219,257,269,272),]");

#Tabulate species by sex
gdf = groupby(penguins, [:species, :sex])
tt = combine(gdf, nrow);
show(IOContext(stdout, :limit=>false), MIME"text/plain"(), tt)

X=penguins[:,3:6]       # select the numerical columns
rics=R"ICS($X, S1=ICS_tcov, S2=ICS_cov)"
ics = ICSModel(S1=tcov, S2=cov2);
scores = fit_predict!(ics, X);


@test(isapprox(ics.kurtosis_, rcopy(rics)[Symbol("gen_kurtosis")]))
@test(isapprox(ics.skewness_, rcopy(rics)[Symbol("gen_skewness")]))
@test(isapprox(ics.W_, rcopy(rics)[Symbol("W")]))
@test(isapprox(scores, rcopy(rics)[Symbol("scores")]))

scree_plot(ics)

## Convert the categorical arrays species and sex to string arrays
using CategoricalArrays
clusters1 = unwrap.(penguins.species)
clusters2 = unwrap.(penguins.sex)

component_plot2(ics, clusters=clusters1)
component_plot2(ics, clusters=clusters2, select=[4])


using CairoMakie
using PairPlots
select = 1:2
ss = DataFrame(scores, ["IC$i" for i=1:size(scores,2)])

pairplot(scores)

fig=pairplot(ss[clusters1 .== "Adelie", :], 
          ss[clusters1 .== "Chinstrap", :],
          ss[clusters1 .== "Gentoo", :] =>
         (
            PairPlots.Scatter(markersize=10),
            PairPlots.MarginHist()
         ), 
          fullgrid=true)


            
table1 = (;
    x = randn(1000),
    y = randn(1000),
)

table2 = (;
    x = 1 .+ randn(1000),
    y = 2 .+ randn(1000),
    z = randn(1000),
)

fig = pairplot(table1, table2)
                    
fig =pairplot(
    PairPlots.Series(ss[clusters1 .== "Adelie", :], color=Makie.wong_colors(0.5)[1]) => (
        PairPlots.Scatter(markersize=8),
        PairPlots.MarginDensity(
            linewidth=2.5f0
        )
    ),
    PairPlots.Series(ss[clusters1 .== "Chinstrap", :], color=Makie.wong_colors(0.5)[2]) => (
        PairPlots.Scatter(markersize=8),
        PairPlots.MarginDensity(
            linewidth=2.5f0
        )
    ),
    PairPlots.Series(ss[clusters1 .== "Gentoo", :], color=Makie.wong_colors(0.5)[3]) => (
        PairPlots.Scatter(markersize=8),
        PairPlots.MarginDensity(
            linewidth=2.5f0
        )
    ),
)

fig=pairplot(
    ss[clusters1 .== "Adelie", :] => (
        PairPlots.Scatter(markersize=8),
        PairPlots.MarginDensity(),
    ),
    ss[clusters1 .== "Chinstrap", :] => (
        PairPlots.Scatter(markersize=8),
        PairPlots.MarginDensity(),
    ),
    ss[clusters1 .== "Gentoo", :] => (
        PairPlots.Scatter(markersize=8),
        PairPlots.MarginDensity(),
    ),
    fullgrid=true
)
Label(fig[0,:], "Component plot of Pengin data set: tcov-COV", fontsize=18)
fig
save("Penguins_plot3c-julia.png", fig)