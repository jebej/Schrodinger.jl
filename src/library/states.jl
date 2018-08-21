"""
    basis(N, n)

Generate a basis state (a.k.a. Fock or number state) ket \$|n⟩\$, in a Hilbert space of size `N`. Note that the size of the Hilbert space must be *at least* `n+1`. The function `fock` is an alias for `basis`.

Returns a sparse vector.

# Example
```jldoctest
julia> ψ = basis(3,2)
3-d Schrodinger.Ket{SparseVector{Float64,Int64},1} with dimensions 3
1.00∠0°|2⟩
```
"""
function basis(n::Integer,dims::SDims)
    N = prod(dims)
    N>n || throw(ArgumentError("basis level $n is too large for a $N-d space"))
    # Julia is 1-indexed!
    return Ket(SparseVector(N,[n+1],[1.0]),dims)
end

basis(N::Integer, n::Integer) = basis(n::Integer,(N,))

"""
    coherent(N, α, analytic=false)

Generate a coherent state ket \$|α⟩\$, in a Hilbert space of size `N`. To create a coherent density operator, use the `Operator` function: `Operator(coherent(N,n))`.

Two methods can be used for generating a coherent state: via application of a displacment operator on a ground state (the default), or analytically, with the formula

```math
|α⟩ = e^{-\\frac{|α|^2}{2}} \\sum_{n=0}^{N-1} \\frac{α^n}{\\sqrt{n!}} |n⟩.
```

While the operator method will return a normalized ket, the analytic method will not. Both methods converge as `N` gets larger. The analytic method is also much faster, especially for large `N`.

Returns a dense vector.

# Example
```jldoctest
julia> coherent(6,0.4+1im)
6-d Schrodinger.Ket{Array{Complex{Float64},1},1} with dimensions 6
0.60∠68°|1⟩ + 0.56∠0°|0⟩ + 0.46∠136°|2⟩ + 0.29∠-155°|3⟩ + 0.15∠-87°|4⟩ +…
```
"""
function coherent(N::Integer, α::Number, analytic::Bool=false)
    if analytic let a = exp(-abs2(α)/2) # due to julia issue #15276
        x = [a*α^n/sqrtfact(n) for n = 0:N-1]
        return Ket(x,(N,))
    end else
        return Ket(data(displacementop(N,α))[:,1],(N,)) # first column of D(α)
    end
end


"""
    thermal(N, n)

Generate a thermal state density matrix \$ρ_n\$ with particle number `n`, in a Hilbert space of size `N`. A thermal state \$ρ_n\$ is a probabilistic mixture of basis states such that the expectation value of the number operator \$\\hat{n}\$ is `n`. Note that this is true only if \$N≫n\$. The returned density matrix is always normalized.

Returns a sparse matrix.

# Example
```jldoctest
julia> N=5; n=0.2;

julia> ρ = thermal(N,n)
5×5 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with dimensions 5
 0.833441  0.0       0.0        0.0         0.0
 0.0       0.138907  0.0        0.0         0.0
 0.0       0.0       0.0231511  0.0         0.0
 0.0       0.0       0.0        0.00385852  0.0
 0.0       0.0       0.0        0.0         0.000643087

julia> expect(numberop(N),ρ)
0.19935691318327978
```
"""
function thermal(N::Integer, n::Real)
    β = log(1.0/n+1.0)
    rowval = collect(1:N)
    colptr = Vector{Int}(N+1); colptr[1:N] = rowval; colptr[end] = N+1
    nzval  = normalize!([exp(-β*k) for k = 0:N-1],1)
    return Operator(SparseMatrixCSC(N,N,colptr,rowval,nzval),(N,),true)
end


"""
    maxmixed(N)

Generate a maximally mixed density matrix. The maximally mixed state is a mixture of basis states with uniform probability.

Returns a sparse matrix.

# Example
```jldoctest
julia> maxmixed(4)
4×4 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with dimensions 4
 0.25  0.0   0.0   0.0
 0.0   0.25  0.0   0.0
 0.0   0.0   0.25  0.0
 0.0   0.0   0.0   0.25
```
"""
function maxmixed(N::Integer)
    rowval = collect(1:N)
    colptr = Vector{Int}(N+1); colptr[1:N] = rowval; colptr[end] = N+1
    nzval  = Vector{Float64}(N); fill!(nzval, 1/N)
    return Operator(SparseMatrixCSC(N,N,colptr,rowval,nzval),(N,),true)
end


"""
    maxentangled(n,N=2)

Generate a maximally entangled state between `n` `N`-d systems:

```math
|\\phi⟩=\\sum_{j=0}^{N-1}\\frac{1}{\\sqrt{N}}|j⟩^{⊗n}.
```

Tracing out all but one of the entangled systems results in a maximally mixed state.

# Example
```jldoctest
julia> ψ = maxentangled(3,4)
64-d Schrodinger.Ket{SparseVector{Float64,Int64},3} with dimensions 4⊗4⊗4
0.50∠0°|0,0,0⟩ + 0.50∠0°|1,1,1⟩ + 0.50∠0°|2,2,2⟩ + 0.50∠0°|3,3,3⟩

julia> ptrace(ψ,(1,3))
4×4 Schrodinger.Operator{Array{Float64,2},1} with dimensions 4
 0.25  0.0   0.0   0.0
 0.0   0.25  0.0   0.0
 0.0   0.0   0.25  0.0
 0.0   0.0   0.0   0.25
```
"""
function maxentangled(n::Int,N::Int=2)
    c = div(N^n-1,N-1)
    nzind = [1+m*c for m=0:N-1]
    nzval = fill(1/sqrt(N),N)
    return Ket(SparseVector(N^n,nzind,nzval),ntuple(_->N,Val{n}))
end


"""
    ket(state,dims=2)

Generate a state ket from a tuple of basis levels and a tuple of corresponding space dimensions. Note that each space dimension must be larger than the level by *at least* 1. If only an integer is passed to `dims`, all basis levels will have the same dimension.

Returns a sparse vector.

# Example
```jldoctest
julia> ψ = ket((3,0,1),(5,2,3))
30-d Schrodinger.Ket{SparseVector{Float64,Int64},3} with dimensions 5⊗2⊗3
1.00∠0°|3,0,1⟩
```

See also: [`@qb_str`](@ref), for generating qubit states via a bitstring.
"""
ket{N}(state::SDims{N},dim::Int=2) = ket(state,ntuple(_->dim,Val(N)))
function ket{N}(state::SDims{N},dims::SDims{N})
    @inbounds for i = 1:N
        dims[i]>state[i] || throw(ArgumentError("basis level $(state[i]) is too large for a $(dims[i])-d space"))
    end
    n = tensored_sub2ind(dims,state) # 0-based indexing in tindex, returns 1 indexed index
    return Ket(SparseVector(prod(dims),[n],[1.0]),dims)
end


"""
    qb"bitstring"

Generate a qubit state ket from the given bitstring. This string macro can also be invoked via the macro notation `@qb_str "bitstring"`.

Returns a sparse vector.

# Example
```jldoctest
julia> Ψ⁻ = normalize!(qb"01" - qb"10")
4-d Schrodinger.Ket{SparseVector{Float64,Int64},2} with dimensions 2⊗2
0.71∠0°|0,1⟩ + 0.71∠180°|1,0⟩
```

See also: [`ket`](@ref), for generating arbitrary states.
"""
macro qb_str(state::AbstractString)
    quote
        dims = $(ntuple(_->2,length(state)))
        Ket(SparseVector(prod(dims),[$(parse(Int,state,2)+1)],[1.0]),dims)
    end
end
