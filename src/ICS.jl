module ICS

using DataFrames: DataFrame
using Statistics
using LinearAlgebra
using Robustbase
using Printf

import Base: show

greet() = print("Hello World from ICS!")

include("Scatter.jl")
include("ICSModel.jl")
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

       greet

end
