# QuObj Tests
using Base.Test, Schrodinger

@testset "QuObj Basic Math" begin
# Build a few different variable for testing
g = basis(2,0)
e1 = basis(2,1)
ρ = maxmixed(4)
σ = create(4)

@testset "Algebra" begin
    # Test algebra between QuObjects and numbers
    @test +g == g
    @test -e1 == 0-e1
    @test +ρ == ρ
    @test -ρ == 0-ρ
    @test data(g+1) == data(1+g) == [2,1]
    @test data(e1-1) == [-1,0]
    @test data(1-e1) == [1,0]
    @test ρ+1 == 1+ρ == ρ+qeye(4)
    @test ρ-0.25 == 0.25-ρ == qzero(4)
    @test data(2*e1) == [0,2]
end

end
