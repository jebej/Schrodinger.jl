module Schrodinger
using DiffEqBase, OrdinaryDiffEq, Compat

export Operator, Ket, Bra,
    data, dims, isnormalized, dimsmatch, dense, braket,
    ptrace, expect, fidelity, fidelity2, levelprobs, tensor,
    fock, basis, coherent, thermal, maxmixed, maxentangled,
    qzero, qeye, create, destroy, numberop, projectorop,
    displacementop, squeezeop,
    Liouvillian, SchrodingerEvo, LindbladEvo,
    Propagator, SchrodingerProp, LindbladProp,
    sesolve, mesolve, lsolve, psolve, psteady,
    expim, gaussian

include(joinpath("quobj","types.jl"))
include(joinpath("quobj","basicmethods.jl"))
include(joinpath("math","basemath.jl"))
include(joinpath("math","eigen.jl"))
include(joinpath("math","special.jl"))
include(joinpath("math","ptrace.jl"))
include(joinpath("dynamics","liouvillian.jl"))
include(joinpath("dynamics","propagator.jl"))
include(joinpath("dynamics","constructors.jl"))
include(joinpath("dynamics","interface.jl"))
include(joinpath("misc","utils.jl"))
include(joinpath("misc","kron.jl"))
include(joinpath("misc","sparsevec.jl"))
include(joinpath("misc","approxherm.jl"))
include(joinpath("library","operators.jl"))
include(joinpath("library","states.jl"))
include(joinpath("library","drivefuns.jl"))
include(joinpath("library","constants.jl"))

#const CACHE = Dict{UInt64,Matrix{Complex128}}()

# 0.5
#scale!(A::Union{Diagonal,Symmetric,Hermitian},b::Number) = (scale!(A.data,b);A)
end