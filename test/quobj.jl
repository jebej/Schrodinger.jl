# QuObj Tests
using Base.Test, Schrodinger

@testset "QuObj Generation" begin

@testset "Ket Generation" begin
    # Test ket generation from arrays and functions
    @test data(Ket([0,0])) == [0.0,0.0]
    @test dims(Ket([0,1,0])) == (3,)
    @test_throws ArgumentError Ket([0,1],(3,))
    @test_throws ArgumentError basis(2,2)
    @test Ket(sparse([0,1/√2,1im/√2,0]),(2,2)) == 1/√2*(basis(2,0)⊗basis(2,1) + 1im*basis(2,1)⊗basis(2,0))
    @test Ket([0,1,0,0,0,0],(3,2)) == (basis(3,0)⊗basis(2,1))
    α = 1.31+0.11im
    @test coherent(30,α)[6] ≈ exp(-0.5*abs2(α))*α^5/sqrt(factorial(5))
    @test coherent(30,α,false) ≈ coherent(30,α,true)
    @test qb"01" == basis(2,0)⊗basis(2,1)
    @test data(ket((2,1,5,2),(5,3,10,4))).nzind[1] == 1+2+4*5+(4*10)*1+(4*10*3)*2
end

@testset "Bra Generation" begin
    # Test bra generation from arrays and functions
    @test data(Bra([0,0])) == [0.0,0.0]
    @test dims(Bra([0,1,0])) == (3,)
    @test_throws ArgumentError Bra([0,1],(3,))
    @test Bra(sparse([0,1/√2,1im/√2,0]),(2,2)) == 1/√2*(basis(2,0)⊗basis(2,1) - 1im*basis(2,1)⊗basis(2,0))'
    α = 1.31+0.11im
    @test Bra(coherent(30,α))[6] ≈ exp(-0.5*abs2(α))*conj(α)^5/sqrt(factorial(5))
end

@testset "Operator Generation" begin
    # Test operator generation from arrays and functions
    @test data(Operator([0 0;0 0])) == zeros(2,2)
    @test dims(Operator([0 0 0;0 1 0;0 0 0])) == (3,)
    @test_throws ArgumentError Operator([0 0;0 0],(3,))
    @test_throws DimensionMismatch Operator([0 1 0;0 0 0])
    @test (A=maxmixed(4); diag(A) == 0.25*ones(4)&&isdiag(A))
    @test (A=thermal(10,0.45); diag(A)[1:2] ≈ [6.89660888E-1,2.14032689E-1] && isdiag(A))
    @test data(qzero(6,(2,3))) == zeros(6,6)
    @test data(qeye(6,(2,3))) == eye(6)
    @test data(qeye(6,(2,3))) == eye(6)
    @test diag(create(6),-1) == [sqrt(n) for n = 1:5]
    @test destroy(6) == create(6)'
    @test (A=numberop(6); diag(A) == [n for n = 0:5] && isdiag(A))
    @test data(displacementop(3,0.5im)) ≈ [0.88261978 0.43980233im -0.16600070; 0.43980233im 0.64785934 0.62197442im; -0.16600070 0.62197442im 0.76523956]
    @test data(squeezeop(3,0.5im)) ≈ [0.93814834 0 -0.34623359im; 0 1 0; -0.34623359im 0 0.93814834]
    @test (U=rand_unitary(10,(5,2)); U*U' ≈ qeye(10,(5,2)))
end

@testset "Matrix Property Checks" begin
    # Tests matrix property checks
    H = Schrodinger.hermitianize!(rand(ComplexF64,5,5))
    H_s = Schrodinger.hermitianize!(sprand(ComplexF64,5,5,0.1))
    x = sparse([1,2,3,4,5],[3,5,3,1,1],[1+5im,3.4,1im,2.5+2.5im,6.6])
    @test isapproxhermitian(H+full(x)*1E-15)
    @test !isapproxhermitian(H+full(x)*1E-12)
    @test isapproxhermitian(H_s+x*1E-15)
    @test !isapproxhermitian(H_s+x*1E-12)
    U1 = sparse([1,2,4,3],[1,2,3,4],[1,1,1,1])
    U2 = expim(Hermitian(H))
    @test isunitary(U1)
    @test !isunitary(U1+sprand(4,4,0.5)*1E-12)
    @test isapproxunitary(U1)
    @test isapproxunitary(U1+sprand(4,4,0.5)*1E-15)
    @test !isapproxunitary(U1+sprand(4,4,0.5)*1E-12)
    @test isunitary(full(U1))
    @test !isunitary(full(U1+sprand(4,4,0.5)*1E-12))
    @test isapproxunitary(U2)
    @test !isapproxunitary(U2+full(x))
end

end
