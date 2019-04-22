using Compat: dropdims

function ptrace_ref(A::Matrix, out, sysdims::NTuple)
    # Function adapted from a [MATLAB function](http://www.dr-qubit.org/Matlab_code.html) by Toby Cubitt and licensed under GPL2.
    # First, calculate systems, dimensions, etc.
    out = collect(out)
    nsub = length(sysdims) # Number of subsytems
    rdims = reverse(sysdims) # Dimensions of subsystems in reverse order
    keep = setdiff(1:nsub,out) # Indices of subsystems to keep
    dimstrace = prod(sysdims[out])::Int # Total dim of subsystems to trace out
    dimskeep = prod(sysdims[keep])::Int # Total dim of subsystems to keep
    # Reshape density matrix into tensor with one row and one column index
    # for each subsystem, permute traced subsystem indices to the end,
    perm = (nsub + 1) .- [reverse(keep); reverse(keep).-nsub; out; out.-nsub]
    B = permutedims(reshape(A,rdims...,rdims...),perm)
    # Reshape again so that first two indices are row and column
    # multi-indices for kept subsystems and third index is a flattened index
    # for traced subsystems
    C = reshape(B,dimskeep,dimskeep,dimstrace^2)
    # Sum third index over "diagonal" entries
    return dropdims(sum(C[:,:,1:dimstrace+1:dimstrace^2],3),dims=3), sysdims[keep]
end

function ptrace_ref(x::Vector, out::NTuple{no,Int}, sysdims::NTuple{ns,Int}) where {no,ns}
    # Make sure the input arguments make sense
    issorted(out) || throw(ArgumentError("subsystem indices $out must be sorted"))
    no>ns && throw(ArgumentError("more subsystem indices ($no) than number of subsystems ($ns)"))
    for i = 1:no
        out[i]>ns && throw(ArgumentError("subsystem index $(out[i]) larger than number of subsystems"))
    end
    # First, calculate systems, dimensions, etc.
    keep = sorted_setdiff(ntuple(identity,Val(ns)),out) # Subsystems to keep
    rsysdims = revtuple(sysdims) # Dimensions of subsystems in reverse order
    keepdims = gettuple(sysdims,keep) # Dimensions of subsystems to keep
    outdims  = gettuple(sysdims,out) # Dimensions of subsystems to trace out
    # Reshape state vector to "reverse" ket on traced subsystems into a bra
    perm = (ns + 1) .- [collect(revtuple(keep)); collect(out)]
    y = permutedims(reshape(x,rsysdims...),perm)
    y = reshape(y,prod(keepdims),prod(outdims))
    # Take outer product
    return y*y', keepdims
end
