using Base: product, tail

"""
    ptrace(ρ, out)

Compute the partial trace of a `Density` matrix or `Operator` ρ by tracing out the subsystems specified by `out`. Multiple subsystems can be traced out by passing a sorted tuple of subsystem indices.

# Example
```jldoctest
Φ₊ = normalize!(basis(2,0)⊗basis(2,0) + basis(2,1)⊗basis(2,1)) # Bell pair
Ψ₊ = normalize!(basis(2,0)⊗basis(2,1) + basis(2,1)⊗basis(2,0)) # Bell pair
ρ  = 0.25 * density(Φ₊) + 0.75 * density(Ψ₊) # density matrix
ptrace(ρ,2) # trace out qubit 2
# output
2×2 Schrodinger.Density{Array{Float64,2},1} with space dimensions 2:
 0.5  0.0
 0.0  0.5
```
"""
function ptrace(ρ::Density, out)
    res = ptrace(full(ρ),out,dims(ρ))
    return Density(res[1],res[2])
end
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
4-d Schrodinger.Ket{SparseVector{Float64,Int64},2} with space dimensions 2⊗2:
0.71∠0°|0,0⟩ + 0.71∠0°|1,1⟩

julia> ptrace(Φ₊,1) # trace out qubit 1
2×2 Schrodinger.Density{Array{Float64,2},1} with space dimensions 2:
 0.5  0.0
 0.0  0.5
```
"""
function ptrace(ψ::Ket, out)
    res = ptrace(full(ψ),out,dims(ψ))
    return Density(res[1], res[2])
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
    # Calculate dimensions and indices
    keep  = sorted_setdiff(ntuple(identity,Val{ns}),out) # Indices of subsystems to keep
    rkeep = revinds(keep,ns)
    rout  = revinds(out,ns)
    rsysdims  = revtuple(sysdims) # Dimensions of subsystems in reverse order
    rkeepdims = gettuple(rsysdims,rkeep)
    routdims  = gettuple(rsysdims,rout)
    # Generate tuples of indices to loop over
    R = product(range.(1,rkeepdims)...)
    S = product(range.(1,routdims)...)
    # Initialize two vectors for indexing purpose
    rinds_ii = Vector{Int}(ns); rinds_jj = Vector{Int}(ns)
    # Initialize output matrix
    B = zeros(T,prod(rkeepdims),prod(rkeepdims))
    # Main loop
    @inbounds for r in R, q in R
        # For each element in the reduced matrix
        i, j = tindexr(q,rkeepdims), tindexr(r,rkeepdims)
        for k=1:(ns-no)
            rinds_ii[rkeep[k]] = q[k]; rinds_jj[rkeep[k]] = r[k]
        end
        for s in S
            # Sum over all the corresponding "diagonal" elements in the full density matrix
            for o=1:no
                rinds_ii[rout[o]] = s[o]; rinds_jj[rout[o]] = s[o]
            end
            ii, jj = tindexr(rinds_ii,rsysdims), tindexr(rinds_jj,rsysdims)
            B[i,j] += _ptrace_ii_jj(A,ii,jj)
        end
    end
    return B, gettuple(sysdims,keep)
end

# Dense partial trace kernel
@inline _ptrace_ii_jj(A::Matrix,ii,jj) = A[ii,jj]
@inline _ptrace_ii_jj(x::Vector,ii,jj) = x[ii]*conj(x[jj])


# ptrace supporting functions
function tindexr{N}(rinds,rsysdims::NTuple{N,Int})
    i = rinds[1]
    d = 1
    @inbounds for n = 2:N
        d *= rsysdims[n-1]
        i += d * (rinds[n]-1)
    end
    return i
end

@inline function sorted_setdiff(t1::Tuple, t2::Tuple)
    if t1[1] == t2[1]
        sorted_setdiff(tail(t1), tail(t2))
    else
        (t1[1], sorted_setdiff(tail(t1), t2)...)
    end
end
@noinline sorted_setdiff(t1::Tuple{}, t2::Tuple) = throw(ArgumentError("duplicate or missing index $(t2[1])"))
sorted_setdiff(t1::Tuple, ::Tuple{}) = t1
sorted_setdiff(::Tuple{}, ::Tuple{}) = ()

revtuple{N}(t::NTuple{N,Any}) = ntuple(i->t[N+1-i],Val{N})
revinds{N}(t::NTuple{N,Any},ns::Int) = ntuple(i->ns+1-t[N+1-i],Val{N})
gettuple{N}(t1::NTuple,t2::NTuple{N,Any}) = ntuple(i->t1[t2[i]],Val{N})
