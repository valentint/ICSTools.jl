using ICS
using Test

@testset "ICS.jl" begin
    @testset "Scatters" begin
    using Robustbase
    X = hbk[:,1:3]
    
    cc=cov2(X)
    @test isapprox(cc.location, [3.2066666666666666, 5.597333333333332, 7.230666666666666])
    @test isapprox(cc.scatter, [13.341711711711707 28.46920720720721 41.243981981981975; 28.46920720720721 67.88296576576586 94.66562342342341; 41.243981981981975 94.66562342342341 137.8348576576576])

    
    end
end

