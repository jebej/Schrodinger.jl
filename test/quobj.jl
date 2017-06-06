# QuObj Tests
using Base.Test, Schrodinger

@testset "QuObj Generation" begin

@testset "Ket Generation" begin
    # Test manual ket generation (from arrays) and library functions
    @test data(Ket([0,0])) == [0.0,0.0]
    @test dims(Ket([0,1,0])) == (3,)
    @test Ket([0,1,0,0],(2,2)) == (basis(2,0)⊗basis(2,1))
    α = 1.31+0.11im
    @test coherent(30,α)[6] ≈ exp(-0.5*abs2(α))*α^5/sqrt(factorial(5))
end

end
