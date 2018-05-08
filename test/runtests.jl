# Schrodinger.jl Tests
using Base.Test, Schrodinger

@testset "All Tests" begin
include("quobj.jl")
include("math.jl")
include("special.jl")
include("dynamics.jl")
include("misc.jl")
end
