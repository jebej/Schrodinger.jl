using Base.product

"""
    ptrace(ρ, out)

Compute the partial trace of a linear `Operator` ρ by tracing out the subsystems specified by `out`. Multiple subsystems can be traced out by passing a sorted tuple of subsystem indices.

# Example
```jldoctest
Φ₊ = normalize!(basis(2,0)⊗basis(2,0) + basis(2,1)⊗basis(2,1)) # Bell pair
Ψ₊ = normalize!(basis(2,0)⊗basis(2,1) + basis(2,1)⊗basis(2,0)) # Bell pair
ρ  = 0.25 * Operator(Φ₊) + 0.75 * Operator(Ψ₊) # density matrix
ptrace(ρ,2) # trace out qubit 2
# output
2×2 Schrodinger.Operator{Array{Float64,2},1} with dimensions 2
 0.5  0.0
 0.0  0.5
```
"""
function ptrace(σ::Operator, out)
    res = ptrace(full(σ),out,dims(σ))
    return Operator(res[1],res[2])
end


"""
    ptrace(ψ, out)

Compute the partial trace of a state `Ket` or `Bra` ψ by tracing out the subsystems specified by `out`. Returns a density matrix. Multiple subsystems can be traced out by passing a sorted tuple of subsystem indices.

# Example
```jldoctest
julia> Φ₊ = normalize!(basis(2,0)⊗basis(2,0) + basis(2,1)⊗basis(2,1)) # Bell pair
4-d Schrodinger.Ket{SparseVector{Float64,Int64},2} with dimensions 2⊗2
0.71∠0°|0,0⟩ + 0.71∠0°|1,1⟩

julia> ptrace(Φ₊,1) # trace out qubit 1
2×2 Schrodinger.Operator{Array{Float64,2},1} with dimensions 2
 0.5  0.0
 0.0  0.5
```
"""
function ptrace(ψ::Ket, out)
    res = ptrace(full(ψ),out,dims(ψ))
    return Operator(res[1], res[2])
end
ptrace(ψ::Bra, out) = ptrace(ψ',out)


# Method for when a single index is passed
ptrace(x, out::Int, sysdims::NTuple) = ptrace(x,(out,),sysdims)

# Dense partial trace
function ptrace{T,no,ns}(A::AbstractArray{T}, out::NTuple{no,Int}, sysdims::NTuple{ns,Int})
    # Make sure the input arguments make sense
    issorted(out) || throw(ArgumentError("subsystem indices $out must be sorted"))
    no>ns && throw(ArgumentError("more subsystem indices ($no) than number of subsystems ($ns)"))
    for i = 1:no
        out[i]>ns && throw(ArgumentError("subsystem index $(out[i]) larger than number of subsystems"))
    end
    # Calculate dimensions
    keep  = sorted_setdiff(ntuple(identity,Val{ns}),out)
    keepdims = gettuple(sysdims,keep)
    outdims  = gettuple(sysdims,out)
    # Generate tuples of subscripts to loop over
    R = product((Base.OneTo.(keepdims).-1)...)
    S = product((Base.OneTo.(outdims).-1)...)
    # Initialize vectors for indexing purpose & output matrix
    subs_ii = Vector{Int}(ns); subs_jj = Vector{Int}(ns)
    B = zeros(T,prod(keepdims),prod(keepdims))
    # Main loop
    @inbounds for r in R, q in R
        # For each element in the reduced matrix
        i, j = tensored_sub2ind(keepdims,q), tensored_sub2ind(keepdims,r)
        for k=1:(ns-no)
            subs_ii[keep[k]] = q[k]; subs_jj[keep[k]] = r[k]
        end
        for s in S
            # Sum over all the corresponding "diagonal" elements in the full density matrix
            for o=1:no
                subs_ii[out[o]] = s[o]; subs_jj[out[o]] = s[o]
            end
            ii, jj = tensored_sub2ind(sysdims,subs_ii), tensored_sub2ind(sysdims,subs_jj)
            B[i,j] += _ptrace_ii_jj(A,ii,jj)
        end
    end
    return B, keepdims
end

# Dense partial trace kernel
@inline _ptrace_ii_jj(A::Matrix,ii,jj) = A[ii,jj]
@inline _ptrace_ii_jj(x::Vector,ii,jj) = x[ii]*x[jj]'
