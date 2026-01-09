using ICSTOOL
using Test
using Random

@testset "ICSTOOL.jl" begin
    @testset "Scatters" begin
    using Robustbase
    X = hbk[:,1:3]
    
    cc=cov2(X);
    display(cc)
    show(cc)
    @test isapprox(cc.location, [3.2066666666666666, 5.597333333333332, 7.230666666666666])
    @test isapprox(cc.scatter, [13.341711711711707 28.46920720720721 41.243981981981975; 28.46920720720721 67.88296576576586 94.66562342342341; 41.243981981981975 94.66562342342341 137.8348576576576])

    cc=covW(X);
    @test isapprox(cc.location, [3.2066666666666666, 5.597333333333332, 7.230666666666666])
    @test isapprox(cc.scatter, [89.73060413073529 238.52394554488217 293.1656684595206; 238.52394554488217 698.7372357715304 805.351453455683; 293.1656684595206 805.351453455683 990.106804274464])
    
    cc=covAxis(X);
    @test isapprox(cc.location, [3.2066666666666666, 5.597333333333332, 7.230666666666666])
    @test isapprox(cc.scatter, [12.652049320295678 25.83111047215079 38.308361936614055; 25.83111047215079 60.28709739852082 86.66570889600766; 38.308361936614055 86.66570889600766 129.2221167907351])

    cc=cov4(X);
    @test isapprox(cc.location, [3.2066666666666666, 5.597333333333332, 7.230666666666666])
    @test isapprox(cc.scatter, [17.94612082614706 47.704789108976435 58.63313369190412; 47.704789108976435 139.74744715430612 161.0702906911366; 58.63313369190412 161.07029069113662 198.02136085489283])

    cc=tcov(X);
    @test isapprox(cc.scatter, [0.6748433227658608 0.5832916311399278 0.9230795838785872; 0.5832916311399278 2.103393485424369 1.9557718887400475; 0.9230795838785872 1.9557718887400475 3.402408498939923])

    Random.seed!(1234)
    cc=mcd_raw(X);
    @test isapprox(cc.location, [1.5333333333333334, 2.456410256410257, 1.6076923076923075])
    @test isapprox(cc.scatter, [2.8170040282000657 0.01100879025915085 0.5147128728712407; 0.01100879025915085 0.8908332251214733 0.43140397241955347; 0.5147128728712407 0.43140397241955347 2.156540524350635])
    
    Random.seed!(1234)
    cc=mcd_raw(X, nsamp="deterministic");
    @test isapprox(cc.location, [1.61025641025641, 2.394871794871795, 1.6743589743589742])
    @test isapprox(cc.scatter, [2.9397144914660727 -0.06904061351203296 1.0144863859568412; -0.06904061351203296 0.9134579667861435 0.19534610988153564; 1.0144863859568412 0.19534610988153564 2.1432468908301505])
    
    let err = nothing
        try
            cc=mcd_raw(X, nsamp="xxx");
        catch err
        end
        @test err isa Exception
        @test sprint(showerror, err) == "Invalid 'nsamp': can be either 'deterministic' or a number!"
    end

    Random.seed!(1234)
    cc=mcd_rwt(X);
    @test isapprox(cc.location, [1.558333333333333, 1.8033333333333335, 1.6599999999999997])
    @test isapprox(cc.scatter, [1.21312099789689 0.02391541790703198 0.1657932538217424; 0.02391541790703198 1.228356794706198 0.19573474773133615; 0.1657932538217424 0.19573474773133615 1.1253468443530048])

    Random.seed!(1234)
    cc=mcd_rwt(X, nsamp="deterministic");
    @test isapprox(cc.location, [1.5377049180327866, 1.780327868852459, 1.6868852459016392])
    @test isapprox(cc.scatter, [1.2208968424415276 0.05473721495128016 0.1265444496168565; 0.05473721495128016 1.242702149861826 0.15178261962426998; 0.1265444496168565 0.15178261962426998 1.1541431351037628])
    
    let err = nothing
        try
            cc=mcd_rwt(X, nsamp="xxx");
        catch err
        end
        @test err isa Exception
        @test sprint(showerror, err) == "Invalid 'nsamp': can be either 'deterministic' or a number!"
    end

    cc=mlc(X, alg=:alg1);
    @test isapprox(cc.location, [1.6089942229177603, 1.8875566256231056, 1.697445667934874])
    @test isapprox(cc.scatter, [1.0790561468851574 0.525861733182914 0.8861365315136714; 0.525861733182914 1.9860137071985495 1.8536943735943225; 0.8861365315136714 1.8536943735943225 3.2277904870231016])

    cc=mlc(X, alg=:alg2);
    @test isapprox(cc.location, [1.6089944684180306, 1.8875565779548642, 1.6974456466045007])
    @test isapprox(cc.scatter, [1.0790542144871826 0.5258609121407963 0.8861350894035047; 0.5258609121407963 1.986010286964204 1.8536914124840267; 0.8861350894035047 1.8536914124840267 3.227785018207513])

    cc=mlc(X);  # alg=:alg3
    @test isapprox(cc.location, [1.6089944684180306, 1.8875565779548642, 1.6974456466045007])
    @test isapprox(cc.scatter, [1.0790542144871826 0.5258609121407963 0.8861350894035047; 0.5258609121407963 1.986010286964204 1.8536914124840267; 0.8861350894035047 1.8536914124840267 3.227785018207513])

    let err = nothing
        try
            cc=mlc(X, alg=:ALG4);
        catch err
        end
        @test err isa Exception
        @test sprint(showerror, err) == "Unknown algorithm: ALG4"
    end

    end
    
    @testset "ICSModel" begin
    using Robustbase
    X = hbk[:,1:3]

    ics = ICSModel();               # by default algorithm="whiten"
    display(ics);                   # model is not fitted yet
    let err = nothing               # model is not fitted yet
        try
            scree_plot(ics);
        catch err
        end
        @test err isa Exception
        @test sprint(showerror, err) == "Model is not fitted yet!"
    end
    let err = nothing               # model is not fitted yet
        try
            outlier_plot(ics);
        catch err
        end
        @test err isa Exception
        @test sprint(showerror, err) == "Model is not fitted yet!"
    end
    let err = nothing               # model is not fitted yet
        try
            component_plot2(ics);
        catch err
        end
        @test err isa Exception
        @test sprint(showerror, err) == "Model is not fitted yet!"
    end

    scores = fit_predict!(ics, X);
    display(ics)
    @test isapprox(ics.kurtosis_, [23.465141777305032, 4.708104428408411, 2.917680272932387])
    @test isapprox(ics.skewness_, [0.01805498369427433, 0.33956728286246396, 0.012932251292713692])

    # all plots
    clusters = repeat(["normal"], outer=size(scores,1))
    clusters[1:14] = repeat(["outlier"], outer=14)
    scree_plot(ics)
    scree_plot(ics, type="bars")
    scree_plot(ics, type="line")
    let err = nothing               # model is not fitted yet
        try
            scree_plot(ics, type="abcd")
        catch err
        end
        @test err isa Exception
        @test sprint(showerror, err) == "Undefined type of plot: abcd!"
    end
    outlier_plot(ics, clusters=clusters)
    component_plot2(ics)
    let err = nothing               # wrong components selected
        try
            component_plot2(ics, select=[1, 2, 3])
        catch err
        end
        @test err isa Exception
        @test sprint(showerror, err) == "Invalid columns selected: should be one or two!"
    end
    let err = nothing               # wrong components selected
        try
            component_plot2(ics, select=[1, 5])
        catch err
        end
        @test err isa Exception
        @test sprint(showerror, err) == "Invalid columns selected: both should be greater than 1 and less than 3"
    end
    let err = nothing               # identical components selected
        try
            component_plot2(ics, select=[2, 2])
        catch err
        end
        @test err isa Exception
        @test sprint(showerror, err) == "Identical columns selected!"
    end

    let err = nothing               # identical components selected
        try
            component_plot2(ics, select=[5])
        catch err
        end
        @test err isa Exception
        @test sprint(showerror, err) == "Invalid column selected: should be greater than 1 and less than 3"
    end

    ics = ICSModel(algorithm="standard");
    scores = fit_predict!(ics, X);
    @test isapprox(ics.kurtosis_, [23.465141777305032, 4.708104428408411, 2.917680272932387])
    @test isapprox(ics.skewness_, [0.01805498369427433, 0.33956728286246396, 0.012932251292713692])

    ics = ICSModel(algorithm="QR");
    scores = fit_predict!(ics, X);
    @test isapprox(ics.kurtosis_, [23.465141777305032, 4.708104428408411, 2.917680272932387])
    @test isapprox(ics.skewness_, [0.01805498369427433, 0.33956728286246396, 0.012932251292713692])

    ics = ICSModel(algorithm="QR", fix_signs="W");
    scores = fit_predict!(ics, X);
    @test isapprox(ics.kurtosis_, [23.465141777305032, 4.708104428408411, 2.917680272932387])
    ## there is no skewness, if fix_signs=="W"
    ##  @test isapprox(ics.skewness_, [0.01805498369427433, 0.33956728286246396, 0.012932251292713692])

    Random.seed!(1234)
    X = randn(300, 9)
    ics = ICSModel()
    fit_predict!(ics, X)
    display(ics)

end
end