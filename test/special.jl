# QuObj Tests
using Base.Test, Schrodinger

@testset "Quantum Information Function Test" begin
# Build a few different variable for testing
N = 13
α = complex(1.1,-0.33)
g = basis(2,0)
e1 = basis(2,1)
numop = numberop(N)

@testset "Expectation Value Tests" begin
    # Test calculation of expectation values
    @test expect(σz,g) == 1
    @test expect(e1,σz) == -1
    @test expect(σx,normalize!(g+e1)) ≈ 1
    @test expect(normalize!(g-e1),σx) ≈ -1
    @test expect(normalize!(g+1im*e1),σy) ≈ 1
    @test expect(σy,normalize!(g-1im*e1)) ≈ -1
    @test expect(coherent(N,α),destroy(N)) ≈ α
    @test expect(destroy(N),coherent(N,0.5α,true)) ≈ 0.5α
    @test [expect(basis(N,i),numop) for i=0:N-1] == collect(0:N-1)
    @test expect(normalize!(basis(N,2)+basis(N,3)),numop) ≈ 2.5
end

end
