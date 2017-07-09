# QuObj Tests
using Base.Test, Schrodinger

@testset "QuObj Basic Math" begin
# Build a few different variable for testing
g = basis(2,0)
e1 = dense(basis(2,1))
ρ = maxmixed(4)
σ = dense(create(4))

@testset "Algebra" begin
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
    @test (x=ρ/0;isnan(x[2,1])&&isinf(x[1,1]))
    @test (x=σ/0;isnan(x[1,1])&&isinf(x[2,1]))
end

end
