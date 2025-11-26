function scree_plot(model::ICSModel; type::String="dots", figsize=(400, 400))

    if isnothing(model.kurtosis_)
        error("Model is not fitted yet!")
    end
    n = size(model.kurtosis_,1)

    ## Create a tidy DataFrame for StatsPlots
    df = DataFrame(
        ID=["IC$i" for i=1:n], 
        Kurtosis=model.kurtosis_,
        index = 1:n
    )

    if(type == "bars")
        bar(df.ID, df.Kurtosis; legend=false, orientation=:v, color=:gray, size=figsize)
    elseif type == "line"
        @df df plot(:index, :Kurtosis, 
            seriestype=:line,
            color=:black, 
            ylab="Generalized kurtosis", 
            legend=false,
            size=figsize
        )
    elseif type == "dots"
    @df df plot(:ID, :Kurtosis, 
        seriestype=:scatter,
        color=:black, 
        ylab="Generalized kurtosis", 
        legend=false,
        size=figsize
    )
    else
        error("Undefined type of plot: ", type, "!")
    end
    
    return current()
end

function outlier_plot(model::ICSModel; clusters::Union{Vector{String}, Nothing}=nothing, legend::Symbol=:topright, figsize=(400, 400))

    if isnothing(model.kurtosis_)
        error("Model is not fitted yet!")
    end
    n = size(model.scores_,1)
    clusters = if isnothing(clusters) repeat(["normal"], outer=n) else clusters end
    ngroup = length(unique(clusters))

    ## Create a tidy DataFrame for StatsPlots
    df = DataFrame(
        ID=1:size(model.scores_, 1), 
        Type=clusters, 
        Z=model.scores_[:,1] .^2
    )

    @df df plot(:ID, :Z, 
        seriestype=:scatter,
        groups=:Type,
        mc=[:lightblue :orange], 
        markershape=[:circle :utriangle], 
        lab=["normal" "outlier"],
        xlab="Observation number", ylab="ICSQR-ICSD2", legendtitle="Type",
        legend=if ngroup <= 1 false else legend end,
        legend_column=-1,
        size=figsize
    ) 
    return current()
end

function component_plot2(model::ICSModel; select=[1, 2], clusters::Union{Vector{String}, Nothing}=nothing, legend::Symbol=:topright, figsize=(400, 400))

    if isnothing(model.kurtosis_)
        error("Model is not fitted yet!")
    end
    n = size(model.scores_, 1)
    p = size(model.scores_, 2)
    names = ["IC$i" for i=1:p]
    clusters = if isnothing(clusters) repeat([1], outer=n) else clusters end
    ngroup = length(unique(clusters))

    if(length(select) == 1)
        if(select[1] < 1 || select[1] > p)
            error("Invalid column selected: should be greater than 1 and less than ", p)
        end
    elseif(length(select) == 2) 
        if(select[1] < 1 || select[1] > p || select[2] < 1 || select[2] > p)
            error("Invalid columns selected: both should be greater than 1 and less than ", p)
        end
        if(select[1] == select[2])
            error("Identical columns selected!")
        end
    else
        error("Invalid columns selected: should be one or two!")
    end
    
    if(length(select) == 1)

        ## Create a tidy DataFrame for StatsPlots
        df = DataFrame(
            X1 = model.scores_[:, select[1]],
            index    = 1:n,
            clusters = clusters
        )
        xlab = "Observation index"
        ylab = names[select[1]]

        @df df plot(:index, :X1,
            seriestype=:scatter,
            groups=:clusters,
    #        mc=[:lightblue :orange], 
    #        markershape=[:circle :utriangle], 
            xlab=xlab, ylab=ylab,
            legend=if ngroup > 1 legend else false end,
            legend_column=-1,
            size=figsize
        ) 
    else

        ## Create a tidy DataFrame for StatsPlots
        df = DataFrame(
            X1 = model.scores_[:, select[1]],
            X2 = model.scores_[:, select[2]],
            index    = 1:n,
            clusters = clusters
        )
        xlab = names[select[1]]
        ylab = names[select[2]]

        @df df plot(:X1, :X2, 
            seriestype=:scatter,
            groups=:clusters,
    #        mc=[:lightblue :orange], 
    #        markershape=[:circle :utriangle], 
            xlab=xlab, ylab=ylab,
            legend=if ngroup > 1 legend else false end,
            legend_column=-1,
            size=figsize
        ) 
    end
    
    return current()
end