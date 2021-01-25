# Schrodinger.jl Tests
using Schrodinger
using Test, LinearAlgebra, SparseArrays

@testset "All Tests" begin
include("quobj.jl")
include("math.jl")
include("special.jl")
include("dynamics.jl")
include("superops.jl")
include("control.jl")
include("tomography.jl")
include("misc.jl")
end
