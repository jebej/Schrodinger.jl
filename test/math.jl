# QuObj Tests
using Base.Test, Schrodinger

@testset "QuObj Basic Math" begin
# Build a few different variable for testing
g = basis(2,0)
e1 = dense(basis(2,1))
ρ = maxmixed(4)
σ = dense(create(4))
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
    @test_broken (x=(g/0);isnan(x[2])&&isinf(x[1])) # see julia PR #22715
    @test (x=(e1/0);isnan(x[1])&&isinf(x[2]))
    @test ρ/2 == qeye(4)/8
    if VERSION>=v"0.6.0"
        @test (x=ρ/0;isnan(x[2,1])&&isinf(x[1,1]))
    else
        @test_broken (x=ρ/0;isnan(x[2,1])&&isinf(x[1,1])) # remove once 0.5 is dropped
    end
    @test (x=σ/0;isnan(x[1,1])&&isinf(x[2,1]))
    @test_throws ArgumentError 3/ρ
    @test_throws ArgumentError 3/σ
    @test ρ^3 == 0.015625*qeye(4)
    @test σx^2 == σy^2 == σz^2 == σ0
    @test_throws MethodError data(ρ^2.5) == data(ρ)^2.5
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
    @test_throws ArgumentError g*e1
    @test ρ*ρ == ρ^2
    @test σ*σ == σ^2
    @test σ*destroy(4) ≈ numberop(4)
    @test destroy(4)'*destroy(4) ≈ numberop(4)
    @test create(4)*σ' ≈ numberop(4)
    @test destroy(4)'*σ' ≈ numberop(4)
    @test σx*g == e1
    @test create(2)'*e1 == g
    @test (g')*destroy(2) == e1'
    @test (g')*create(2)' == e1'
    @test g'*g == (e1')*e1 == 1
    @test dot(plus,plus) ≈ 1
end

end
