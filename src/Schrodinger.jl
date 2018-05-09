module Schrodinger
using Base: tail, front, promote_eltype
using DiffEqBase, OrdinaryDiffEq, Optim, Compat

export Operator, Ket, Bra,
    data, dims, isnormalized, dimsmatch, dense, braket, isapproxhermitian,
    isunitary, isapproxunitary,
    ptrace, expect, fidelity, fidelity2, levelprobs, tensor,
    fock, basis, coherent, thermal, maxentangled, maxmixed, ket, qb,
    qzero, qeye, create, destroy, numberop, projectorop,
    displacementop, squeezeop,
    rand_unitary,
    Liouvillian, SchrodingerEvo, LindbladEvo,
    Propagator, SchrodingerProp, LindbladProp,
    sesolve, mesolve, lsolve, psolve, psteady,
    grape, NormPSU, CoherentSubspaces,
    expim, gaussian, inner

include("quobj/types.jl")
include("quobj/basicmethods.jl")
include("math/basemath.jl")
include("math/eigen.jl")
include("math/special.jl")
include("math/ptrace.jl")
include("dynamics/liouvillian.jl")
include("dynamics/propagator.jl")
include("dynamics/constructors.jl")
include("dynamics/interface.jl")
include("control/grape.jl")
include("control/normpsu.jl")
include("control/coherentsubspaces.jl")
include("control/amplitudelimiter.jl")
include("misc/utils.jl")
include("misc/inner.jl")
include("misc/kron.jl")
include("misc/sparsevec.jl")
include("misc/checks.jl")
include("misc/hermfact.jl")
include("library/operators.jl")
include("library/states.jl")
include("library/random.jl")
include("library/drivefuns.jl")
include("library/constants.jl")

# A few VERSION-conditional definitions
include("basepatch/v0.6.jl")

end
