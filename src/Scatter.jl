"""
A struct to represent the scatter matrix and its related data.

# Fields:
    location::Union{Vector{Float64}, Nothing}: The mean location of the data.
    scatter::Matrix{Float64}: The scatter matrix.
    label::String: A label describing the scatter matrix.
"""
struct Scatter
    location::Union{Vector{Float64}, Nothing}
    scatter::Matrix{Float64}
    label::String
end

function Base.show(io::IO, obj::Scatter)
    print_object(io, obj, multiline = false)
end

function Base.show(io::IO, mime::MIME"text/plain", obj::Scatter)
    #   you can add IO options if you want
    #multiline = get(io, :multiline, true)
    #print_object(io, obj, multiline = multiline)

    println(io, "-> Scatter: " , obj.label)
    
    if !isnothing(obj.location)
        println("Location:")
        println(IOContext(stdout, :compact=>true), obj.location)
    end

    println()
    println(io, "Scatter:")
    Base.show(stdout, mime, obj.scatter)
end

function print_object(io::IO, obj::Scatter; multiline::Bool)
    if multiline
        show(io, MIME"text/plain", obj)
    else
        # write something short, or go back to default mode
        Base.show_default(io, obj)
    end
end

"""
    cov2(X::Union{Matrix{Float64}, DataFrame}; location::Bool=true)

Compute the covariance matrix.

# Arguments:
    X::Union{Matrix{Float64}, DataFrame}: The data matrix.
    location::Bool (default=true): Whether to include the mean location.

# Returns:
    Scatter: An object containing the location and scatter matrix.
"""
function cov2(X::Union{Matrix{Float64}, DataFrame}; location::Bool=true)
    
    if X isa DataFrame
        X = Matrix(X)
    end

    scatter_ = cov(X, dims=1)  # covariance along columns
    location_ = location ? vec(mean(X, dims=1)) : nothing
    return Scatter(location_, scatter_, "COV")
end

"""
    covW(X::Union{Matrix{Float64}, DataFrame}; location::Bool=true, alpha=1, cf=1)

Estimates the scatter matrix based on one-step M-estimator.

# Arguments:
    X::Union{Matrix{Float64}, DataFrame}: The data matrix.
    location::Bool (default=true): Whether to include the mean location.
    alpha (default=1): Parameter of the one-step M-estimator.
    cf (default=1): Consistency factor of the one-step M-estimator.

# Returns:
    Scatter: An object containing the location and weighted scatter matrix.
"""
function covW(X::Union{Matrix{Float64}, DataFrame}; location::Bool=true, alpha=1, cf=1)
    if X isa DataFrame
        X = Matrix(X)
    end

    n, p = size(X)

    if any(ismissing.(X))
        throw(ArgumentError("Missing values are not allowed in X"))
    end
    if p <= 1
        throw(ArgumentError("X must be at least bi-variate"))
    end

    ## Calculate mean and covariance
    X_means = vec(mean(X, dims=1))
    X_cov = cov(X, dims=1)

    ## Mahalanobis distance squared
    inv_cov = inv(X_cov)
    dists = [dot(x .- X_means, inv_cov * (x .- X_means)) for x in eachrow(X)]
    weights = dists .^ alpha

    ## Center data
    X_centered = X .- X_means'
    X_covW = cf / n * (X_centered' * Diagonal(weights) * X_centered)

    location_ = location ? X_means : nothing
    return Scatter(location_, X_covW, "COVW")
end

"""
    covAxis(X::Union{Matrix{Float64}, DataFrame}; location::Bool=true)

Compute the one-step Tyler shape matrix (CovAxis).

# Arguments:
    X::Union{Matrix{Float64}, DataFrame}: The data matrix.
    location::Bool (default=true): Whether to include the mean location.

# Returns:
    Scatter: An object containing the location and scatter matrix.
"""
function covAxis(X::Union{Matrix{Float64}, DataFrame}; location::Bool=true)
    if X isa DataFrame
        X = Matrix(X)
    end

    p = size(X, 2)
    covaxis_scatter = covW(X; location=location, alpha=-1, cf=p)
    return Scatter(covaxis_scatter.location, covaxis_scatter.scatter, "COVAxis")
end

"""
    cov4(X::Union{Matrix{Float64}, DataFrame}; location::Bool=true)

Compute a custom weighted covariance matrix (cov4).

# Arguments:
    X::Union{Matrix{Float64}, DataFrame}: The data matrix.
    location::Bool (default=true): Whether to include the mean location.

# Returns:
    Scatter: An object containing the location and custom weighted scatter matrix.
"""
function cov4(X::Union{Matrix{Float64}, DataFrame}; location::Bool=true)
    if X isa DataFrame
        X = Matrix(X)
    end

    p = size(X, 2)
    cov4_scatter = covW(X; location=location, alpha=1, cf=1 / (p + 2))
    return Scatter(cov4_scatter.location, cov4_scatter.scatter, "COV4")
end


"""
    tcov(X::Union{Matrix{Float64}, DataFrame}, beta=2)

Computes a pairwise one-step M-estimate of scatter with weights based on pairwise 
Mahalanobis distances. 

Note that it is based on pairwise differences and therefore 
does not require a location estimate.

# Arguments:
X::Union{Matrix{Float64}, DataFrame}: The data matrix.
beta: A positive numeric value specifying the tuning parameter of the pairwise one-step M-estimator (defaults to 2)


# Returns:
Scatter: An object containing the location= and pairwise one-step M-estimate of scatter.

# References
Caussinus, H. and Ruiz-Gazen, A. (1993) Projection Pursuit and Generalized Principal Component Analysis. In Morgenthaler, S., Ronchetti, E., Stahel, W.A. (eds.) New Directions in Statistical Data Analysis and Robustness, 35-46. Monte Verita, Proceedings of the Centro Stefano Franciscini Ascona Series. Springer-Verlag.

Caussinus, H. and Ruiz-Gazen, A. (1995) Metrics for Finding Typical Structures by Means of Principal Component Analysis. In Data Science and its Applications, 177-192. Academic Press.

-------------
Compute the TCOV scatter matrix (pairwise, exponentially weighted) without parallelization.

This is a direct Julia translation/optimization of the Armadillo code:
- Uses a Cholesky factor of the sample covariance to *whiten* data for fast Mahalanobis distances.
- Accumulates the scatter in the **original** space using BLAS symmetric rank-1 updates (`syr!`).
- Falls back to an inverse-based path if the covariance is not positive definite.

The weight is `w = exp((-beta/2) * r^2)`, where `r^2` is the squared Mahalanobis distance of a pair.
"""
function tcov(X::Union{Matrix{Float64}, DataFrame}, beta=2)
    if X isa DataFrame
        X = Matrix(X)
    end

    n, p = size(X)
    n < 2 && throw(ArgumentError("X must have at least two rows"))

    # Sample covariance of columns
    Σ = cov(X)
    b = -float(beta) / 2.0

    # Preallocate accumulators
    V = zeros(p, p)
    denom = 0.0
    d_orig = zeros(p)

    # Try Cholesky path for fast Mahalanobis distances
    try
        F = cholesky(Σ) # Σ = U'U
        # Whitened data: Xw = X * inv(U)
        Xw = X / F.U # triangular solve, no explicit inverse
        d_w = similar(d_orig)

        @inbounds for i in 2:n
            xi = @view X[i, :]
            xiw = @view Xw[i, :]
            for j in 1:(i-1)
                xj = @view X[j, :]
                xjw = @view Xw[j, :]

                # Differences in original and whitened spaces
                d_orig .= xi .- xj
                d_w .= xiw .- xjw

                # Squared Mahalanobis distance via whitened diff
                r2 = sum(abs2, d_w)
                w = exp(b * r2)

                # Symmetric rank-1 update: V += w * d_orig * d_orig'
                BLAS.syr!('U', w, d_orig, V)
                denom += w
            end
        end
    catch e
        if e isa PosDefException
            # Fallback: explicit inverse (slower, but robust)
            Σ⁻¹ = inv(Σ)
            @inbounds for i in 2:n
                xi = @view X[i, :]
                for j in 1:(i-1)
                    xj = @view X[j, :]
                    d_orig .= xi .- xj
                    r2 = dot(d_orig, Σ⁻¹ * d_orig)
                    w = exp(b * r2)
                    BLAS.syr!('U', w, d_orig, V)
                    denom += w
                end
            end 
        else
            rethrow()
        end
    end

    # Symmetrize accumulated upper triangle and normalize
    return Scatter(nothing, Matrix(Symmetric(V, :U)) / denom, "TCOV")

end

"""
    mcd_raw(X::Union{Matrix{Float64}, DataFrame}; 
        location::Bool=true, alpha=0.5, nsamp::Union{Int, String})

Compute a raw MCD estimate.

# Arguments:
    X::Union{Matrix{Float64}, DataFrame}: The data matrix.
    location::Bool (default=true): Whether to include the mean location.
    alpha: 
    nsamp:

# Returns:
    Scatter: An object containing the location and custom weighted scatter matrix.
"""
function mcd_raw(X::Union{Matrix{Float64}, DataFrame}; 
        location::Bool=true, alpha=0.5, nsamp::Union{Int, String}=500)

    if nsamp isa String && nsamp == "deterministic"
        mcd = Robustbase.DetMcd(alpha=alpha, reweighting=false)
        Robustbase.fit!(mcd, X)
    elseif nsamp isa Number
        mcd = Robustbase.CovMcd(alpha=alpha, n_initial_subsets=nsamp, reweighting=false)
        Robustbase.fit!(mcd, X)
    else
        error("Invalid 'nsamp': can be either 'deterministic' or a number!")        
    end    
    
    if location 
        location = Robustbase.location(mcd) 
    else 
        location = nothing
    end
    
    return Scatter(location, Robustbase.covariance(mcd), "MCD")
end

"""
    mcd_rwt(X::Union{Matrix{Float64}, DataFrame}; 
        location::Bool=true, alpha=0.5, nsamp::Union{Int, String})

Compute a reweighted MCD estimate.

# Arguments:
    X::Union{Matrix{Float64}, DataFrame}: The data matrix.
    location::Bool (default=true): Whether to include the mean location.
    alpha: 
    nsamp:

# Returns:
    Scatter: An object containing the location and custom weighted scatter matrix.
"""
function mcd_rwt(X::Union{Matrix{Float64}, DataFrame}; 
        location::Bool=true, alpha=0.5, nsamp::Union{Int, String}=500)

    if nsamp isa String && nsamp == "deterministic"
        mcd = Robustbase.DetMcd(alpha=alpha, reweighting=true)
        Robustbase.fit!(mcd, X)
    elseif nsamp isa Number
        mcd = Robustbase.CovMcd(alpha=alpha, n_initial_subsets=nsamp, reweighting=true)
        Robustbase.fit!(mcd, X)
    else
        error("Invalid 'nsamp': can be either 'deterministic' or a number!")        
    end    

    if location 
        location = Robustbase.location(mcd) 
    else 
        location = nothing
    end
    
    return Scatter(location, Robustbase.covariance(mcd), "RMCD")
end
