# QuObj Tests
using Schrodinger
using Schrodinger: ComplexF64
using Compat.Test, Compat.LinearAlgebra, Compat.SparseArrays
println("Testing Math...")

@testset "QuObj Basic Math" begin
# Build a few different variable for testing
g = basis(2,0)
e1 = dense(basis(2,1))
ρ = maxmixed(4)
α = rand(ComplexF64)
ψα = coherent(15,α)
β = rand(ComplexF64)
ψβ = coherent(15,β)
σ = dense(create(4))
a4 = destroy(4)
a2 = destroy(2)
adag2 = create(2)
plus = normalize!(g+1.0im*basis(2,1))

@testset "QuObj/Number Algebra" begin
    # Test algebra between QuObjects and numbers
    @test +g == g
    @test -g == 0-g
    @test +e1 == e1
    @test -e1 == 0-e1
    @test +(g') == g'
    @test -(g') == 0-g'
    @test +(e1') == e1'
    @test -(e1') == 0-e1'
    @test +ρ == ρ
    @test -ρ == 0-ρ
    @test +σ == σ
    @test -σ == 0-σ
    @test g+1 == 1+g == Ket([2,1])
    @test e1-1 == Ket([-1,0])
    @test 1-e1 == Ket([1,0])
    @test g'+1 == 1+g' == Bra([2,1])
    @test e1'-1 == Bra([-1,0])
    @test 1-e1' == Bra([1,0])
    @test ρ+1 == 1+ρ == ρ+qeye(4)
    @test ρ-0.25 == 0.25-ρ == qzero(4)
    @test σ+2 == 2+σ == σ+2*qeye(4)
    @test σ-2 == -(2-σ) == σ-2*qeye(4)
    @test 2*g == g*2 == Ket(sparse([2,0]))
    @test 2*e1 == e1*2 == Ket([0,2])
    @test 2*(g') == (g')*2 == Bra(sparse([2,0]))
    @test 2*(e1') == (e1')*2 == Bra([0,2])
    @test 4*ρ == ρ*4 == qeye(4)
    @test 2σ == σ*2 == σ+σ
    @test g/2 == Ket(sparse([0.5,0]))
    @test e1/1 == e1
    @test (g')/2 == Bra(sparse([0.5,0]))
    @test (e1')/1 == e1'
    @test_throws ArgumentError 2/g
    @test_throws ArgumentError 2/e1
    @test_throws ArgumentError 2/(g')
    @test_throws ArgumentError 2/(e1')
    if VERSION > v"0.7.0-"
        @test (x=(g/0);isnan(x[2])&&isinf(x[1])) # see julia PR #22715
    end
    @test (x=(e1/0);isnan(x[1])&&isinf(x[2]))
    @test ρ/2 == qeye(4)/8
    @test (x=ρ/0;isnan(x[2,1])&&isinf(x[1,1]))
    @test (x=σ/0;isnan(x[1,1])&&isinf(x[2,1]))
    @test_throws ArgumentError 3/ρ
    @test_throws ArgumentError 3/σ
    @test_throws ArgumentError plus^3
    @test_throws ArgumentError plus^3.3
    @test ρ^3 == 0.015625*qeye(4)
    @test σx^2 == σy^2 == σz^2 == σ0
    if VERSION < v"0.7.0-"
        @test_throws MethodError data(ρ^2.5)
    else
        @test data(ρ^2.5) == data(ρ)^2.5
    end
    @test data(dense(ρ)^2.5) == full(ρ)^2.5
    @test (σ^2)[4,2] == √(2)*√(3)
    @test_broken data(σ^2.5) == data(σ)^2.5
end

@testset "QuObj/QuObj Algebra" begin
    @test g+g == 2g
    @test g-g == Ket(sparse([0,0]))
    @test e1+e1+e1 == 3*e1
    @test g+e1 == Ket([1,1])
    @test g-e1 == Ket([1,-1])
    @test g'+g' == 2g'
    @test e1'+e1'+e1' == 3*e1'
    @test g'+e1' == Bra([1,1])
    @test g'-e1' == Bra([1,-1])
    @test ρ+ρ == 2ρ
    @test σ+σ+σ == 3σ
    @test σ-σ == qzero(4)
    @test ρ-ρ == qzero(4)
    @test data(ρ+σ) == data(ρ)+data(σ)
    @test data(ρ-σ) == data(ρ)-data(σ)
    @test_throws ArgumentError g+σ
    @test_throws ArgumentError g-σ
    @test ρ*ρ == ρ^2
    @test σ*σ == σ^2
    @test σ*a4 ≈ numberop(4)
    @test a4'*a4 ≈ numberop(4)
    @test σ*σ' ≈ numberop(4)
    @test a4'*σ' ≈ numberop(4)
    @test adag2*g == adag2*Bra(g)'== e1
    @test adag2'*e1 == adag2'*Bra(e1)' == g
    @test g'*a2 == Bra(g)*a2 == e1'
    @test g'*adag2' == Bra(g)*adag2' == e1'
    @test g'*g == Bra(e1)*e1 == 1
    @test dot(plus,plus) ≈ 1
    @test dot(ψα,ψβ) ≈ exp(-(abs2(α)+abs2(β))/2+α'*β)
    @test plus*plus' == plus*Bra(plus) ≈ normalize!(σ0 + σy)
    @test_throws ArgumentError g*e1
    @test_throws ArgumentError g/e1
end

@testset "Operator Inner Product" begin
    @test inner(σx,σx) == 2
    @test inner(σx,σy) == 0
    @test inner(σx,σz) == 0
    @test inner(ρ*projectorop(size(ρ,1),1:2),ρ) == 2/size(ρ,1)^2
    @test inner(Operator(ψα),Operator(ψβ)) ≈ exp(-abs2(α-β))
    @test_throws DimensionMismatch inner(Operator(ψα),ρ)
end

@testset "Tensor Product" begin
    @test g ⊗ g   == Ket(Bra(g)  ⊗ Bra(g)) == Ket(sparse([1,0,0,0]),(2,2))
    @test g ⊗ e1  == Ket(Bra(g)  ⊗ Bra(e1)) == Ket(sparse([0,1,0,0]),(2,2))
    @test e1 ⊗ g  == Ket(Bra(e1) ⊗ Bra(g)) == Ket(sparse([0,0,1,0]),(2,2))
    @test e1 ⊗ e1 == Ket(Bra(e1) ⊗ Bra(e1)) == Ket([0,0,0,1],(2,2))
    @test data(Operator(g) ⊗ σ) == [data(σ) zeros(4,4); zeros(4,8)]
    @test data(ρ ⊗ Operator(g)) == (A=Diagonal(fill(0.25,8));A.diag[2:2:8].=0;A)
    @test_throws ArgumentError e1 ⊗ ρ
    @test_throws ArgumentError Bra(e1) ⊗ σ
end

@testset "Transposition" begin
    @test g' == Bra([1,0])
    @test (g+1im*e1)' == Bra([1,-1im])
    @test conj((g+1im*e1)) == Ket([1,-1im])
    @test Bra(g)' == Ket([1,0])
    @test Bra([1,1im])' == Ket([1,-1im])
    @test conj(Bra([1,1im])) == Bra([1,-1im])
    @test ρ' == transpose(ρ) == ρ
    @test σ' == transpose(σ)
end

@testset "Math Functions" begin
    for f in [sqrt,log,exp], A in [Operator(coherent(4,1.2+3im)),ρ,create(4)+destroy(4)]
        @test data(f(A)) == f(full(A))
    end
end

@testset "Misc Functions" begin
    for f in [real,imag,abs,abs2], A in [g,e1,plus,Bra(plus),σ,ρ,Operator(coherent(4,1.2+3im))]
        @test f(A) == f.(data(A))
    end
end

@testset "Eigen" begin
    E = 3.3; Δ = 0.25
    O = dense(E*σz+Δ*σx)
    @test eigvals(O) == [-√(E^2+Δ^2), √(E^2+Δ^2)]
    @test eigen(O)[2][1] ≈ normalize!(Ket([-(E-√(E^2+Δ^2))/Δ,-1]))
    a = destroy(20) ⊗ qeye(2)
    s = qeye(20) ⊗ σ₋
    fc = 5; fa = 6; g = 0.025;
    H = fc*(a'a+0.5) + fa*(s's-0.5) + g*(a's + a*s')
    function En(l)
        n,q,Δ = (l-1)÷2, 1-2*(l%2), fa-fc
        return l==0 ? -Δ/2 : (n+1)*fc + 0.5q*√(4g^2*(n+1)+Δ^2)
    end
    θ(n) = 0.5*atan(2g*√(n+1)/(fa-fc))
    D,V = eigs(H,nev=10)
    @test D ≈ En.(0:9)
    @test abs(inner(V[6], cos(θ(2))*ket((3,0),(20,2)) - sin(θ(2))*ket((2,1),(20,2)))) ≈ 1
end

end
