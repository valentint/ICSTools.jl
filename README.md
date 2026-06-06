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
The package ```ICSTools``` can be installed using the ```Julia``` ```REPL``` as follows

```
julia> ]
(@v1.12) pkg> add ICSTools
```

or

```
julia> import Pkg
julia> Pkg.add("ICSTools")
```

and then

```
julia> using ICSTools
```
to make it available to the user.

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

# Citation 
If you use this software, please cite:

Todorov, V. (2026). Tools for Exploring Multivariate Data: The Invariant Coordinate Selection (ICS) Method.
URL https://valentint.github.io/ICSTools.jl/dev/. doi:10.5281/zenodo.20476499.


or in ```BibTeX``` format:

```
@software{ICSTools.jl,
    author = {Valentin Todorov},
    doi = {10.5281/zenodo.20476499},
    license = {["MIT"]},
    title = {{Tools for Exploring Multivariate Data: The Invariant Coordinate Selection (ICS) Method}},
    year = 2026,
    url = {https://valentint.github.io/ICSTools.jl/dev/}}
```

## Community guidelines

### Report issues and request features

If you experience any bugs or issues or if you have any suggestions for
additional features, please submit an issue via the
[*Issues*](https://github.com/valentint/ICSTools.jl/issues) tab of this
repository. Please have a look at existing issues first to see if your
problem or feature request has already been discussed.

### Contribute to the package

If you want to contribute to the package, you can fork this repository
and create a pull request after implementing the desired functionality.

### Ask for help

If you need help using the package, or if you are interested in
collaborations related to this project, please get in touch with the
package maintainer.
