module ICS

using DataFrames: DataFrame
using Statistics
using LinearAlgebra
using Robustbase
using Printf
using StatsPlots

import Base: show

include("Scatter.jl")
include("ICSModel.jl")
include("plots.jl")
include("utils.jl")

export Scatter,
       cov2,
       covW,
       cov4,
       covAxis,
       tcov,
       mcd_raw,
       mcd_rwt,
       
       ICSModel,
       fit!,
       predict,
       fit_predict!,
       scree_plot,
       outlier_plot

end
