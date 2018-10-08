module Schrodinger
using Compat
using Compat.LinearAlgebra, Compat.SparseArrays, Compat.Printf
using Base: tail, front, product, promote_eltype
using OrdinaryDiffEq, Optim
using OrdinaryDiffEq.DiffEqBase: AbstractParameterizedFunction

export Operator, Ket, Bra,
    data, dims, isnormalized, dimsmatch, dense, braket, isapproxhermitian,
    isunitary, isapproxunitary, tensor, ptrace, expect, levelprobs,
    fidelity, fidelity2, entanglement_fidelity, gate_fidelity,
    fock, basis, coherent, thermal, maxentangled, maxmixed, ket, @qb_str,
    qzero, qeye, create, destroy, numberop, projectorop,
    displacementop, squeezeop,
    rand_unitary,
    Liouvillian, SchrodingerEvo, LindbladEvo,
    Propagator, SchrodingerProp, LindbladProp,
    sesolve, mesolve, lsolve, psolve, psteady,
    grape, NormPSU, CoherentSubspaces,
    expim, gaussian, inner,
    eigen

# VERSION-conditional definitions
if VERSION < v"1.1.0-"
    include("basepatch/v1.0.jl")
end
if VERSION < v"0.7.0-"
    include("basepatch/v0.6.jl")
end
if VERSION > v"0.7.0-"
    import Base: BitSet
    import LinearAlgebra: adjoint, exp, exp!, tr, eigen
    import Arpack: eigs
    export eigs, trace
    trace(A::AbstractMatrix) = tr(A)
end

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
include("misc/tupletools.jl")
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

end
