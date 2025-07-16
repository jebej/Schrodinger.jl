module Schrodinger
using LinearAlgebra, SparseArrays, Printf
using Base: tail, front, product
using OrdinaryDiffEq, Optim

# Export submodules
export Gate

# Export types & functions
export Operator, Ket, Bra,
    data, dims, isnormalized, dimsmatch, dense, braket, isapproxhermitian,
    isunitary, isapproxunitary, isdensityop, checkdensityop,
    hermitian, symmetric,
    tensor, ptrace, expect, levelprobs, findstate,
    purity, fidelity, fidelity2, entanglement_fidelity, gate_fidelity,
    fock, basis, coherent, thermal, maxentangled, maxmixed, ket, @qb_str,
    qzero, qeye, create, destroy, numberop, projectorop,
    displacementop, squeezeop,
    rand_unitary,
    Liouvillian, SchrodingerEvo, LindbladEvo,
    Propagator, SchrodingerProp, LindbladProp,
    sesolve, mesolve, lsolve, psolve, psteady,
    grape, NormPSU, CoherentSubspaces,
    expim, gaussian, inner, scale, scale!

# re-export some LinearAlgebra methods
export eigen, eigs, eigvals, trace, normalize, normalize!

include("basepatch/v1.0.jl")
@static if VERSION < v"1.1"
    using Future: copy!
end

include("quobj/types.jl")
include("quobj/basicmethods.jl")
include("math/basemath.jl")
include("math/eigen.jl")
include("math/special.jl")
include("math/superops.jl")
include("math/norms.jl")
include("math/ptrace.jl")
include("dynamics/liouvillian.jl")
include("dynamics/propagator.jl")
include("dynamics/constructors.jl")
include("dynamics/interface.jl")
include("control/grape.jl")
include("control/normpsu.jl")
include("control/coherentsubspaces.jl")
include("control/amplitudelimiter.jl")
include("tomography/statetomo.jl")
include("tomography/processtomo.jl")
include("misc/utils.jl")
include("misc/indexing.jl")
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
include("library/gates.jl")

include("quobj/gpu.jl")

end
