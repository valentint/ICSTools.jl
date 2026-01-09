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

"""
    mlc(X::Union{Matrix{Float64}, DataFrame}; location::Bool=true, 
        alg::Symbol=:alg3, mu_init=nothing, V_init=nothing, gamma_init=nothing, 
        eps=1e-6, maxiter=100)


    Compute Cauchy location and scatter estimates
    It is a wrapper for the Cauchy estimator of location and scatter for a 
        multivariate t-distribution, as computed by tM().

# Arguments:
    X::Union{Matrix{Float64}, DataFrame}: The data matrix.
    location::Bool (default=true): Whether to include the mean location.
    alg::Symbol (default=:alg3): specifies which algorithm to use. Options are :alg1, :alg2 or :alg3. 
    mu_init: initial value for the location vector if available.
    V_init: initial value for the scatter matrix if available.
    gamma_init: initial value for gamma if available. Only needed for alg2.
    eps: convergence tollerance.
    maxiter: maximum number of iterations.

# Returns:
    Scatter: An object containing the location (if reqested) and a numeric 
        matrix giving the estimate of the scatter matrix.

# Details:
This function implements the EM algorithms described in Kent et al. (1994). 
    The norm used to define convergence is as in Arslan et al. (1995).

Algorithm 1 is valid for all degrees of freedom df > 0. 
Algorithm 2 is well defined only for degrees of freedom df > 1. 
Algorithm 3 is the limiting case of Algorithm 2 with degrees of freedom df = 1.

The performance of the algorithms are compared in Arslan et al. (1995).

# References:
Kent, J.T., Tyler, D.E. and Vardi, Y. (1994), A curious likelihood identity for the multivariate t-distribution, Communications in Statistics, Simulation and Computation, 23, 441–453. <doi:10.1080/03610919408813180>.

Arslan, O., Constable, P.D.L. and Kent, J.T. (1995), Convergence behaviour of the EM algorithm for the multivariate t-distribution, Communications in Statistics, Theory and Methods, 24, 2981–3000. <doi:10.1080/03610929508831664>.

"""
function mlc(X::Union{Matrix{Float64}, DataFrame}; location::Bool=true, 
    alg::Symbol=:alg3, mu_init=nothing, V_init=nothing, gamma_init=nothing, eps=1e-6, maxiter=100)

    if X isa DataFrame
        X = Matrix(X)
    end

    n, p = size(X)

    mu_init === nothing && (mu_init = vec(mean(X, dims=1)))
    V_init === nothing && (V_init = cov(X))

    ## we fix the df to have only cauchy estimate
    mu, V = tM(X; df=1, alg=alg, mu_init=mu_init, V_init=V_init, 
        gamma_init=gamma_init, eps=eps, maxiter=maxiter)

    location_ = location ? mu : nothing
    return Scatter(location_, V, "MLC")

end

"""
    tM(X; df::Real=1, alg::Symbol=:alg3, mu_init=nothing, V_init=nothing, 
        gamma_init=nothing, eps=1e-6, maxiter=100)


Compute joint M-estimation of Location and Scatter for a Multivariate t-distribution.
    Implements three EM algorithms to M-estimate the location vector and scatter matrix 
    of a multivariate t-distribution.

# Arguments:
    X::Union{Matrix{Float64}, DataFrame}: The data matrix.
    df::Real (default=1): assumed degrees of freedom of the t-distribution. 
        Default is 1 which corresponds to the Cauchy distribution.
    alg::Symbol (default=:alg3): specifies which algorithm to use. Options are :alg1, :alg2 or :alg3. 
    mu_init: initial value for the location vector if available.
    V_init: initial value for the scatter matrix if available.
    gamma_init: initial value for gamma if available. Only needed for alg2.
    eps: convergence tollerance.
    maxiter: maximum number of iterations.

# Returns:
    mu: vector with the estimated loaction.
    V:  matrix of the estimated scatter.
    gam: estimated value of gamma. Only present when alg2 is used.
    iter: number of iterations.

# Details:
This function implements the EM algorithms described in Kent et al. (1994). 
    The norm used to define convergence is as in Arslan et al. (1995).

Algorithm 1 is valid for all degrees of freedom df > 0. 
Algorithm 2 is well defined only for degrees of freedom df > 1. 
Algorithm 3 is the limiting case of Algorithm 2 with degrees of freedom df = 1.

The performance of the algorithms are compared in Arslan et al. (1995).

# References:
Kent, J.T., Tyler, D.E. and Vardi, Y. (1994), A curious likelihood identity for the multivariate t-distribution, Communications in Statistics, Simulation and Computation, 23, 441–453. <doi:10.1080/03610919408813180>.

Arslan, O., Constable, P.D.L. and Kent, J.T. (1995), Convergence behaviour of the EM algorithm for the multivariate t-distribution, Communications in Statistics, Theory and Methods, 24, 2981–3000. <doi:10.1080/03610929508831664>.

"""
function tM(X; df::Real=1, alg::Symbol=:alg3, 
    mu_init=nothing, V_init=nothing, gamma_init=nothing, eps=1e-6, maxiter=100)

    n, p = size(X)

    mu_init === nothing && (mu_init = vec(mean(X, dims=1)))
    V_init === nothing && (V_init = cov(X))

    if alg != :alg2 && gamma_init !== nothing
        @warn "Initial gamma is only used in alg2"
    end

    gamma_init === nothing && (gamma_init = 1.0)

    if alg == :alg1
        return alg1(X, mu_init, V_init, df, eps, maxiter)
    elseif alg == :alg2
        return alg2(X, mu_init, V_init, gamma_init, df, eps, maxiter)
    elseif alg == :alg3
        return alg3(X, mu_init, V_init, df, eps, maxiter)
    else
        error("Unknown algorithm: $alg")
    end
end


function norm_mu_V(a::AbstractVector, B::AbstractMatrix, A::AbstractMatrix)
    Ainv = inv(A)
    BAinv = B * Ainv
    square =
        dot(a, Ainv * a) +
        sum(diag(BAinv * BAinv))
    return sqrt(square)
end

function alg1(X, mu_init, V_init, nu, eps, maxiter)
    n, p = size(X)

    V = copy(V_init)
    mu = copy(mu_init)

    iter = 0
    differ = Inf

    while differ > eps
        iter += 1

        d2 = Robustbase.mahalanobis_distance(X, mu, V) .^ 2
        u = (nu + p) ./ (nu .+ d2)

        mu_new = vec(mean(X .* u, dims=1)) / mean(u)

        Xc = X .- mu_new'
        V_new = (Xc' * (Xc .* u)) / n

        differ = norm_mu_V(mu_new - mu, V_new - V, V_new)

        mu = mu_new
        V = V_new

        if iter ≥ maxiter
            error("maxiter reached without convergence")
        end
    end

    return (mu=mu, V=V, iter=iter)
end

function alg2(X, mu_init, V_init, gamma_init, nu, eps, maxiter)
    n, p = size(X)

    V = copy(V_init)
    mu = copy(mu_init)
    gamma = gamma_init

    iter = 0
    differ = Inf

    while differ > eps
        iter += 1

        d2 = Robustbase.mahalanobis_distance(X, mu, V) .^ 2
        w = (nu + p) ./ (nu - 1 + 1/gamma .+ d2 ./ gamma)

        gamma_new = mean(w)
        mu_new = vec(mean(X .* w, dims=1)) / gamma_new

        Xc = X .- mu_new'
        V_new = (Xc' * (Xc .* w)) / n / gamma_new

        differ = norm_mu_V(mu_new - mu, V_new - V, V_new)

        gamma = gamma_new
        mu = mu_new
        V = V_new

        if iter ≥ maxiter
            error("maxiter reached without convergence")
        end
    end

    return (mu=mu, V=V, gam=gamma, iter=iter)
end

function alg3(X, mu_init, V_init, nu, eps, maxiter)
    n, p = size(X)

    V = copy(V_init)
    mu = copy(mu_init)
    gamma = 1.0

    iter = 0
    differ = Inf

    while differ > eps
        iter += 1

        d2 = Robustbase.mahalanobis_distance(X, mu, V) .^ 2
        w = (nu + p) ./ (nu - 1 + 1/gamma .+ d2 ./ gamma)

        gamma_new = mean(w)
        mu_new = vec(mean(X .* w, dims=1)) / gamma_new

        Xc = X .- mu_new'
        V_new = (Xc' * (Xc .* w)) / n / gamma_new

        differ = norm_mu_V(mu_new - mu, V_new - V, V_new)

        mu = mu_new
        V = V_new
        gamma = gamma_new

        if iter ≥ maxiter
            error("maxiter reached without convergence")
        end
    end

    return (mu=mu, V=V, iter=iter)
end

