## Assume Scatter, cov, covW, covAxis, cov4, 
##  sort_eigenvalues_eigenvectors, sqrt_symmetric_matrix, 
##  sign_max, check_gen_kurtosis are already defined elsewhere.

mutable struct ICSModel
    S1::Function
    S2::Function
    algorithm::String
    center::Bool
    fix_signs::String
    S1_args::Dict{Symbol,Any}
    S2_args::Dict{Symbol,Any}
    W_::Union{Nothing,Matrix{Float64}}
    scores_::Union{Nothing,Matrix{Float64}}
    kurtosis_::Union{Nothing,Vector{Float64}}
    skewness_::Union{Nothing,Vector{Float64}}
    feature_names_in_::Union{Nothing,Vector{String}}
    S1_X_::Union{Nothing,Any}  # will hold a Scatter
end

"""
    ICSModel(; S1=cov, S2=covW, algorithm="whiten", center=false,
              fix_signs="scores", S1_args=Dict(), S2_args=Dict())

Constructor for an ICS model.

Invariant Coordinate Selection (ICS) Class and associated methods.

The ICSModel class implements the ICS algorithm: it transforms the data, 
via the simultaneous diagonalization of two scatter matrices, into 
an invariant coordinate system or independent components, depending 
on the underlying assumptions.

It supports various scatter matrix calculations and offers multiple 
algorithms for applying ICS.

# Parameters:

    S1 (function returning a scatter object): (default: cov2) Function to compute the first scatter matrix.

    S2 (function returning a scatter object): (default: covW) Function to compute the second scatter matrix.

    algorithm::String: (default: 'whiten') The algorithm used for transformation ('standard', 'whiten', 'QR').

    center::Bool: (default: False): a logical indicating whether the invariant coordinates should be centered with respect to the first locattion or not. Centering is only applicable if the first scatter object contains a location component, otherwise this is set to False. Note that this only affects the scores of the invariant components (attribute scores_), but not the generalized kurtosis values (attribute kurtosis_).

    fix_signs::String: (default: 'scores') How to fix the signs of the invariant coordinates. Possible values are 'scores' to fix the signs based on (generalized) skewness values of the coordinates, or 'W' to fix the signs based on the coefficient matrix of the linear transformation.

    S1_args (dict): Additional arguments for S1.

    S2_args (dict): Additional arguments for S2.

# Returns:
    An object of class ICSModel with the following fields:

    `W_::Matrix{Float64}`: Transformation matrix in which each row contains the coefficients of the linear transformation to the corresponding invariant coordinate.

    `scores_::Matrix{Float64}`: Transformed matrix in which each column contains the scores of the corresponding invariant coordinate.

    `kurtosis_::Vector{Float64}`: Generalized kurtosis values.

    `skewness_::Vector{Float64}`: Skewness values.

    `feature_names_in_::Vector{String}`: Names of features seen during fit. Defined only when X has feature names that are all strings.

    `S1_X_::Any`: Fitted scatter S1. Defined only when center=True.

# Details:

Supported algorithms:

    1. standard: performs the spectral decomposition of the symmetric matrix :math:`S_1(X)^{-1/2}S_2(X)S_1(X)^{-1/2}`

    2. whiten: whitens the data with respect to the first scatter matrix before computing the second scatter matrix.

    3. QR: numerically stable algorithm based on the QR algorithm for a common family of scatter pairs: if S1 is cov(), and if S2 is one of cov4, covW, or covAxis. See Archimbaud et al. (2023) for details.


See also [`cov2`](@ref), [`cov4`](@ref), [`covW`](@ref), [`covAxis`](@ref), [`tcov`](@ref).

# Examples:
    
```julia
using ICSTools, RDatasets

iris = dataset("datasets", "iris")

x = iris[:,1:4]

X = Matrix(iris[:,1:4])

ics = ICSModel()    # default S1=cov2 and S2=covW

fit!(ics, X)        # fit!() modifies its first argument, ics.

```
"""
function ICSModel(; S1=cov2, S2=covW, algorithm="whiten", center=false,
                   fix_signs="scores", S1_args=Dict{Symbol,Any}(),
                   S2_args=Dict{Symbol,Any}())

    valid_algo = ["whiten", "standard", "QR"]
    @assert algorithm in valid_algo "algorithm must be one of $valid_algo"

    valid_fix_signs = ["scores", "W"]
    @assert fix_signs in valid_fix_signs "fix_signs must be one of $valid_fix_signs"

    return ICSModel(S1, S2, algorithm, center, fix_signs,
                    S1_args, S2_args,
                    nothing, nothing, nothing, nothing, nothing, nothing)
end

function Base.show(io::IO, model::ICSModel)
    println(io, "\nICS based on two scatter matrices")
    println(io, "S1: ", nameof(model.S1))
    println(io, "S1_args: ", model.S1_args)
    println(io, "S2: ", nameof(model.S2))
    println(io, "S2_args: ", model.S2_args)

    println(io, "\nInformation on the algorithm:")
    println(io, "algorithm: ", model.algorithm)
    println(io, "center: ", model.center)
    println(io, "fix_signs: ", model.fix_signs)

    if isnothing(model.kurtosis_)
        println()
        println("Model is not fitted yet!")
    else
        p = length(model.kurtosis_)
        if p <= 6
            select = 1:p
        else
            select = vcat(1:3, ((p-2):p))
        end
        println(io, "\nThe generalized kurtosis measures of the components are:")
        if p <= 6
            for (i, val) in enumerate(model.kurtosis_)
                println(io, "IC_$(i): ", @sprintf("%.4f", val))
            end
        else
            for (i, val) in enumerate(model.kurtosis_[1:3])
                println(io, "IC_$(i): ", @sprintf("%.4f", val))
            end
            println(io, "... ", @sprintf("%4s", "    "))
            for (i, val) in enumerate(model.kurtosis_[(p-2):p])
                println(io, "IC_$((p-3)+i): ", @sprintf("%.4f", val))
            end
        end

        println(io, "\nThe coefficient matrix of the linear transformation is:")
        if !isnothing(model.W_)
            feature_names = isnothing(model.feature_names_in_) ?
                            ["Feature_$i" for i in 1:p] :
                            model.feature_names_in_
            if p <= 6
                # header row
                header = "     " * join([@sprintf("%12s", name) for name in feature_names], " ")
                println(io, header)

                for (i, row) in enumerate(eachrow(model.W_))
                    row_str = join([@sprintf("%12.5f", val) for val in row], " ")
                    println(io, @sprintf("IC_%-3d %s", i, row_str))
                end
            else
                # header row
                select = vcat(1:3, (p-2):p)
                header1 = "     " * join([@sprintf("%13s", name) for name in feature_names[1:3]], " ")
                header2 = join([@sprintf("%13s", name) for name in feature_names[(p-2):p]], " ")
                println(io, header1,  @sprintf("%13s", "..."), header2)

                for (i, row) in enumerate(eachrow(model.W_)[1:3])
                    row_str1 = join([@sprintf("%14.5f", val) for val in row[1:3]], " ")
                    row_str2 = join([@sprintf("%14.5f", val) for val in row[(p-2):p]], " ")
                    println(io, @sprintf("IC_%-3d %s %7s %s", i, row_str1, "...", row_str2))
                end
                
                println("...")
                for (i, row) in enumerate(eachrow(model.W_)[(p-2):p])
                    row_str1 = join([@sprintf("%14.5f", val) for val in row[1:3]], " ")
                    row_str2 = join([@sprintf("%14.5f", val) for val in row[(p-2):p]], " ")
                    println(io, @sprintf("IC_%-3d %s %7s %s", p-3+i, row_str1, "...", row_str2))
                end
            end
        else
            println(io, "None")
        end
    end
end

# ------------------- Main methods -------------------
"""
    fit!(model::ICSModel, X::Union{Matrix{Float64}, DataFrame})

Fit the ICS model to data.

This function relies on several helper methods to perform the ICS transformation:
_validate_input, _compute_first_scatter, _compute_second_scatter,
_transform_second_scatter, _compute_transformation, _compute_transformation_qr,
_center_data, _fix_component_signs.

Parameters:

    X::Union{Matrix{Float64}, DataFrame}: Data to fit the ICS model, where rows are samples and columns are features.

Returns:
    
    The fitted ICS object.
"""
function fit!(model::ICSModel, X::Union{Matrix{Float64}, DataFrame})
    if X isa DataFrame
        model.feature_names_in_ = names(X)
        X = Matrix(X)
    end

    n, p = size(X)
    if isnothing(model.feature_names_in_)
        column_names = ["Feature_$i" for i in 1:p]
        model.feature_names_in_ = column_names
    end
    
    _validate_input(X)

    S1_X, S1_X_inv_sqrt = _compute_first_scatter(model, X)

    if model.algorithm == "whiten"
        Y = X * S1_X_inv_sqrt
        S2_Y = _compute_second_scatter(model, Y)
        W, gen_kurtosis = _compute_transformation(model, S1_X_inv_sqrt, S2_Y)
    elseif model.algorithm == "QR"
        W, gen_kurtosis = _compute_transformation_qr(model, X, S1_X)
    else
        S2_X = _compute_second_scatter(model, X)
        S2_Y = _transform_second_scatter(S1_X_inv_sqrt, S2_X)
        W, gen_kurtosis = _compute_transformation(model, S1_X_inv_sqrt, S2_Y)
    end

    W_final, gen_skewness = _fix_component_signs(model, X, W)

    if model.center
        model.S1_X_ = S1_X
    end

    model.W_ = W_final
    model.kurtosis_ = gen_kurtosis
    model.skewness_ = gen_skewness

    return model
end

"""
    predict(model::ICSModel, X::Union{Matrix{Float64}, DataFrame})

Transform the data using the fitted ICS model.

Parameters:

    X::::Union{Matrix{Float64}, DataFrame}: Data to transform.

Returns:

    Matrix{Float64}: Transformed matrix in which each column contains the scores of the corresponding invariant
    coordinate.
"""
function predict(model::ICSModel, X::Union{Matrix{Float64}, DataFrame})
    if X isa DataFrame
        X = Matrix(X)
    end

    isnothing(model.W_) && error("ICS model must be fitted before calling predict().")
    n, p = size(X)
    @assert size(model.W_, 1) == p "Expected $(size(model.W_, 1)) features in X."

    if model.center && !isnothing(model.S1_X_)
        X = _center_data(model, X, model.S1_X_)
    end

    Z_final = X * model.W_'
    model.scores_ = Z_final
    return Z_final
end

"""
    fit_predict!(model::ICSModel, X::Union{Matrix{Float64}, DataFrame}) = predict(fit!(model, X), X)

Fit the ICS model and transform the data using the fitted ICS model.

Parameters:

    X::Union{Matrix{Float64}, DataFrame}: Data to fit and transform.

Returns:

    Matrix{Float64}: Transformed matrix in which each column contains the scores of the corresponding invariant
    coordinate.
"""
fit_predict!(model::ICSModel, X::Union{Matrix{Float64}, DataFrame}) = predict(fit!(model, X), X)

# ------------------- Helpers -------------------

function _validate_input(X::AbstractMatrix)
    n, p = size(X)
    @assert p > 1 "X must be at least bi-variate"
    return true
end

function _compute_first_scatter(model::ICSModel, X)
    S1_X = model.S1(X; model.S1_args...)
    S1_X_inv_sqrt = sqrt_symmetric_matrix(S1_X.scatter; inverse=true)
    return S1_X, S1_X_inv_sqrt
end

function _compute_second_scatter(model::ICSModel, X)
    return model.S2(X; model.S2_args...)
end

function _transform_second_scatter(S1_X_inv_sqrt, S2_X)
    return Scatter(nothing, S1_X_inv_sqrt * S2_X.scatter * S1_X_inv_sqrt, S2_X.label)
end

function _compute_transformation(model::ICSModel, S1_X_inv_sqrt, S2_Y)
    evals, evecs = eigen(Symmetric(S2_Y.scatter))
    evals, evecs = sort_eigenvalues_eigenvectors(evals, evecs)
    W = evecs' * S1_X_inv_sqrt
    return W, evals
end

function _compute_transformation_qr(model::ICSModel, X, S1_X)
    n, p = size(X)
    T1_X = S1_X.location
    X_centered = isnothing(T1_X) ? X : X .- T1_X'

    norms = [maximum(abs, X_centered[i,:]) for i in 1:n]
    order_rows = sortperm(norms, rev=true)
    X_reordered = X_centered[order_rows, :]

    ##  F = qr(X_reordered / sqrt(n-1); pivot=true)
    ##  Q, R, pvt = Matrix(F.Q), F.R, F.prow

    F = qr(X_reordered / sqrt(n-1), ColumnNorm())
    Q, R, pvt = Matrix(F.Q), F.R, F.jpvt

    d = (n-1) .* sum(Q.^2, dims=2)

    if model.S2 == cov4
        α, cf = 1, 1/(p+2)
    elseif model.S2 == covAxis
        α, cf = -1, p
    elseif model.S2 == covW
        α = get(model.S2_args, :alpha, 1)
        cf = get(model.S2_args, :cf, 1)
    end

    d = d[:]
    S2_Y = cf * (n-1)/n * (Q .* (d.^α))' * Q

    evals, evecs = eigen(Symmetric(S2_Y))
    check_gen_kurtosis(evals)
    evals, evecs = sort_eigenvalues_eigenvectors(evals, evecs)

    W = (R \ evecs)'   # solve linear system
    W = W[:, invperm(pvt)]
    return W, evals
end

function _center_data(model::ICSModel, X, S1_X)
    T1_X = S1_X.location
    if isnothing(T1_X)
        @warn "Location component in S1 is required for centering. Proceeding without centering."
        model.center = false
        return X
    else
        return X .- T1_X'
    end
end

function _fix_component_signs(model::ICSModel, X, W)
    if model.fix_signs == "scores"
        Z = X * W'
        gen_skewness = vec(mean(Z, dims=1) .- mapslices(median, Z; dims=1))
        signs = map(x -> x ≥ 0 ? 1.0 : -1.0, gen_skewness)
        gen_skewness = signs .* gen_skewness
        W_final = W .* signs
        return W_final, gen_skewness
    else
        row_signs = [sign_max(W[i,:]) for i in 1:size(W,1)]
        row_norms = sqrt.(sum(W.^2, dims=2))

        ## This is wrong (resulting in differences to R), 
        ##  remove the transpose 
        ##  W_final = (W' ./ (row_signs .* row_norms))'
        
        W_final = (W ./ (row_signs .* row_norms))
        return W_final, nothing
    end
end