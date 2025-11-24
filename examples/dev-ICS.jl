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
##  Change to the root directory of the package, activate the package and run 'test'
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

ICS.outlier_plot(ics_qr, clusters)

df = DataFrame(ID=1:size(scores, 1), Type=clusters, Z=scores[:,1] .^2)

fig = @df df plot(:ID, :Z, 
    seriestype=:scatter,
    groups=:Type,
    mc=[:lightblue :orange], 
    markershape=[:circle :utriangle], 
    lab=["normal" "outlier"],
    xlab="Observation number", ylab="ICSQR-ICSD2", legendtitle="Type" )


savefig(fig, "C:/projects/statproj/julia/HTP3_outlier.png")

##==================================================================
##
##  data_philips from cellWise

using DataFrames
using RCall
using Test

## Load R libraries 
R"library(ICSOutlier)"
R"library(ICSClust)"

R"library(cellWise)"
X = rcopy(R"data('data_philips', package='cellWise'); x=data_philips");

rics=R"ICS($X, S1=ICS_mcd_raw, S2=ICS_cov, S1_args = list(alpha = 0.50, nsamp = 500))"
ics = ICSModel(S1=mcd_raw, S2=cov2);
scores = fit_predict!(ics, X);

@test(isapprox(ics.kurtosis_, rcopy(rics)[Symbol("gen_kurtosis")]))
@test(isapprox(ics.skewness_, rcopy(rics)[Symbol("gen_skewness")]))
@test(isapprox(ics.W_, rcopy(rics)[Symbol("W")]))
@test(isapprox(scores, rcopy(rics)[Symbol("scores")]))

using CairoMakie
using PairPlots

clusters = Vector{String}(undef, size(scores,1))
fill!(clusters, "Group1")
clusters[1:100] .= "Group2"
clusters[491:565] .= "Group3"

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
penguins = X = rcopy(R"data('penguins', package='palmerpenguins'); x=penguins");

using StatsPlots

#=
    using Plots
    plot(
    penguins.bill_length_mm, 
    penguins.bill_depth_mm, 
    seriestype=:scatter,
    size=(500,500)
    )
=#

@df penguins plot(
  :bill_length_mm,
  :bill_depth_mm,
  seriestype=:scatter,
  group=:species,
  size=(500,500)
)

@df penguins plot(
  :bill_length_mm,
  :bill_depth_mm,
  seriestype=:scatter,
  group=:species,
  title="Palmer Penguins",
  xlabel="Bill Length (mm)",
  ylabel="Bill Depth (mm)",
  size=(500,500)
)

using StatsPlots
@df penguins plot(
  :bill_length_mm,
  :bill_depth_mm,
  seriestype=:scatter,
  group=:species,
  markershape=[:circle :diamond :utriangle],
  markercolor=[:red :blue :orange],
  markersize=6,
  title="Palmer Penguins",
  xlabel="Bill Length (mm)",
  ylabel="Bill Depth (mm)",
  size=(500,500)
)