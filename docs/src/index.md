# ICSTools.jl

Documentation for ICSTools.jl

Tools for Exploring Multivariate Data: The Invariant Coordinate Selection (ICS) Method.

# Scatters
```@docs
ICSTools.Scatter
```

## Covariance matrix: cov2
```@docs
ICSTools.cov2
```
## One-step M-estimator: covW
```@docs
ICSTools.covW
```
## One-step Tyler shape matrix: CovAxis
```@docs
ICSTools.covAxis
```

## Custom weighted covariance matrix: cov4
```@docs
ICSTools.cov4
```

## Pairwise one-step M-estimate of scatter: tcov
```@docs
ICSTools.tcov
```

## Raw MCD estimate: mcd_raw
```@docs
ICSTools.mcd_raw
```

## Reweighted MCD estimate: mcd_rwt
```@docs
ICSTools.mcd_rwt
```

## Joint M-estimation of Location and Scatter for a Multivariate t-distribution: tM_base

```@docs
ICSTools.tM_base
```

## Multivariate t location and scatter estimates: tM
```@docs
ICSTools.tM
```

## Cauchy location and scatter estimates: mlc
```@docs
ICSTools.mlc
```

## Aggregated standardised variance/covariance in a Mahalanobis neighbourhood of the data points: lcov
```@docs
ICSTools.lcov
```

# ICS Model

## Constructor for an ICS model
```@docs
ICSTools.ICSModel
```

## Fit the ICS model to data
```@docs
ICSTools.fit!
```

## Transform the data using the fitted ICS model
```@docs
ICSTools.predict
```

## Fit the ICS model and transform the data using the fitted ICS model
```@docs
ICSTools.fit_predict!
```

