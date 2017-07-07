# QuObj Tests
using Base.Test, Schrodinger

@testset "Quantum Information" begin
# Build a few different variable for testing
N = 13
α = complex(1.1,-0.33)
g = basis(2,0)
e1 = basis(2,1)
numop = numberop(N)

@testset "Expectation Value" begin
    # Test calculation of expectation valu
    @test expect(σz,g) == 1
    @test expect(e1,σz) == -1
    @test expect(σx,normalize!(g+e1)) ≈ 1
    @test expect(normalize!(g-e1),σx) ≈ -1
    @test expect(normalize!(g+1im*e1),σy) ≈ 1
    @test expect(σy,normalize!(g-1im*e1)) ≈ -1
    @test expect(coherent(N,α),destroy(N)) ≈ α
    @test expect(destroy(N),coherent(N,0.5α,true)) ≈ 0.5α
    @test [expect(basis(N,i),numop) for i=0:N-1] == collect(0:N-1)
    @test expect(normalize!(basis(N,2)+2*basis(N,3)),numop) == (2+2^2*3)/5
    @test expect(normalize!(Operator(basis(N,2))+2*Operator(basis(N,3))),numop) == (1*2+2*3)/3
    @test expect(numop,thermal(N,0.21)) ≈ 0.21
    @test expect(projectorop(N,[2,4]),maxmixed(N)) == 2/N
end

@testset "Fidelity" begin
    # Test calculation of fidelity
    @test fidelity(g,g) == 1
    @test fidelity2(g,g) == 1
    @test fidelity(g,e1) == 0
    @test fidelity2(g,e1) == 0
    @test fidelity(g,normalize!(g+e1)) == 1/√(2)
    @test fidelity2(g,normalize!(g+e1)) ≈ 1/2
    @test fidelity(Operator(coherent(N,α)),basis(N,3)) ≈ abs(coherent(N,α)[4])
    @test fidelity(basis(N,3),Operator(coherent(N,0.5α))) ≈ abs(coherent(N,0.5α)[4])
    @test fidelity2(Operator(coherent(N,α)),basis(N,3)) ≈ abs2(coherent(N,α)[4])
    @test fidelity2(basis(N,3),Operator(coherent(N,0.5α))) ≈ abs2(coherent(N,0.5α)[4])
    @test fidelity(Operator(coherent(N,α)),Operator(basis(N,3))) ≈ abs(coherent(N,α)[4])
    @test fidelity(Operator(basis(N,2)),Operator(coherent(N,1))) ≈ abs(coherent(N,1)[3])
    @test fidelity2(Operator(coherent(N,α)),Operator(basis(N,3))) ≈ abs2(coherent(N,α)[4])
    @test fidelity2(Operator(basis(N,2)),Operator(coherent(N,1))) ≈ abs2(coherent(N,1)[3])
end

end
