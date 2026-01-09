# ICS.jl

Documentation for ICS.jl

Tools for Exploring Multivariate Data: The InvariantCoordinate Selection (ICS) Method.

# Scatters
```@docs
ICS.Scatter
```

## Covariance matrix: cov2
```@docs
ICS.cov2
```
## One-step M-estimator: covW
```@docs
ICS.covW
```
## One-step Tyler shape matrix: CovAxis
```@docs
ICS.covAxis
```

## Custom weighted covariance matrix: cov4
```@docs
ICS.cov4
```

## Pairwise one-step M-estimate of scatter: tcov
```@docs
ICS.tcov
```

## Raw MCD estimate: mcd_raw
```@docs
ICS.mcd_raw
```

## Reweighted MCD estimate: mcd_rwt
```@docs
ICS.mcd_rwt
```

## Cauchy location and scatter estimates: mlc
```@docs
ICS.mlc
```

## Joint M-estimation of Location and Scatter for a Multivariate t-distribution: tM

```@docs
ICS.tM
```
# ICS Model

## Constructor for an ICS model
```@docs
ICS.ICSModel
```

## Fit the ICS model to data
```@docs
ICS.fit!
```

## Transform the data using the fitted ICS model
```@docs
ICS.predict
```

## Fit the ICS model and transform the data using the fitted ICS model
```@docs
ICS.fit_predict!
```

