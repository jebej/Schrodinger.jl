# Schrodinger.jl Tests
using Base.Test, Schrodinger

@testset "Schrodinger.jl Tests" begin
include("quobj.jl")
include("math.jl")
include("special.jl")
include("misc.jl")
end
