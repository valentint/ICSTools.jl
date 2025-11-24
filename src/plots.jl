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

function outlier_plot(model::ICSModel; clusters::Vector{String}, figsize=(400, 400))

    if isnothing(model.kurtosis_)
        error("Model is not fitted yet!")
    end
    n = size(model.scores_,1)

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
        size=figsize
    ) 
    return current()
end