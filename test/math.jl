# QuObj Tests
using Base.Test, Schrodinger

@testset "QuObj Basic Math" begin
# Build a few different variable for testing
g = basis(2,0)
e1 = basis(2,1)

@testset "Addition" begin
    # Test addition between QuObjects and numbers
    @test data(g+1) == [2,1]
end

end
