# ICSTools

[![Build Status](https://github.com/valentint/ICSTools.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/valentint/ICSTools.jl/actions/workflows/CI.yml?query=branch%3Amain)


[![codecov.io](http://codecov.io/github/valentint/ICSTools.jl/coverage.svg?branch=main)](http://codecov.io/github/valentint/ICSTools.jl?branch=main)
[![Doc](https://img.shields.io/badge/docs-stable-blue.svg)](https://valentint.github.io/ICSTools.jl/stable/)
[![Doc](https://img.shields.io/badge/docs-dev-blue.svg)](https://valentint.github.io/ICSTools.jl/dev/)
[![DOI](https://zenodo.org/badge/1080038478.svg)](https://doi.org/10.5281/zenodo.20476499)

# Overview
Invariant Coordinate Selection (ICS) is a data transformation method using simultaneous diagonalization of two scatter matrices, into an invariant coordinate system or independent components. The transformation depends on the underlying assumptions. It is particularly useful for dimension reduction. Unlike PCA, ICS is not based on variance maximization but on the maximization/minimization of a generalized kurtosis, and it is invariant not only to orthogonal data transformations but to any affine transformation.

The ```ICSTools``` package brings the main functionalities of the ```R``` package ```ICS``` to Julia, offering tools for identifying and selecting invariant coordinates in multivariate data. It provides various scatter estimators, transformation settings, and plotting utilities. Our extensive testing ensures results consistent with the ```R``` package, making it easy for users to transition from ```R``` to ```Julia``` or start fresh with ICS.

Check out the [documentation](https://valentint.github.io/ICSTools.jl/stable/) for more details.

See also the ```Python``` package ```ICSpyLab``` at [GitHub](https://github.com/cbecquart/ICSpyLab?tab=readme-ov-file).

![](Readme-logo.png=250x)<!-- -->

# Installation
```
import Pkg; Pkg.add("ICSTools")
using ICSTools
```

# Usage
```{julia}
using ICSTools

## Load dataset
using Robustbase
X=wood[:,1:5];

## Instantiate ICSModel object with all default parameters
ics = ICSModel()
show(ics)

scores = fit_predict!(ics, X);
show(ics)
scores

## scree plot and 2-dimensional component plot
using Plots
gr()
scree_plot(ics)

component_plot2(ics)                # by default the first two components

component_plot2(ics, select=[3,4])  # select components 3 and 4
```