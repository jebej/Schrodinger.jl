# Schrodinger.jl Tests
using Schrodinger
using Schrodinger: ComplexF64

using Compat.Test, Compat.LinearAlgebra, Compat.SparseArrays

@testset "All Tests" begin
include("quobj.jl")
include("math.jl")
include("special.jl")
include("dynamics.jl")
include("tomography.jl")
include("misc.jl")
end
