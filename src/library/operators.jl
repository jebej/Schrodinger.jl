"""
    qzero(N, dims=(N,))

Generate a zero operator for a Hilbert space of size `N`. It is possible to specify the subspace dimensions with the `dims` argument. Returns a sparse matrix.

# Example
```jldoctest
julia> qzero(4,(2,2))
4×4 Operator{SparseMatrixCSC{Float64,Int64},2} with dimensions 2⊗2
 0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0
```
"""
function qzero(::Type{T},N::Integer, dims::Dims=(N,)) where {T<:Number}
    rowval = Vector{Int}(undef,0)
    colptr = ones(Int,N+1)
    nzval  = Vector{T}(undef,0)
    return Operator(SparseMatrixCSC(N,N,colptr,rowval,nzval),dims,true)
end
qzero(::Type{T},dims::Dims) where {T<:Number} = qzero(T,prod(dims),dims)
qzero(N::Integer, dims::Dims=(N,)) = qzero(Float64,N,dims)
qzero(dims::Dims) = qzero(Float64,dims)

"""
    qeye(N, dims=(N,))

Generate an identity operator for a Hilbert space of size `N`. It is possible to specify the subspace dimensions with the `dims` argument. Returns a sparse matrix.

# Example
```jldoctest
julia> qeye(4,(2,2))
4×4 Operator{SparseMatrixCSC{Float64,Int64},2} with dimensions 2⊗2
 1.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0
 0.0  0.0  1.0  0.0
 0.0  0.0  0.0  1.0
```
"""
function qeye(N::Integer, dims::Dims=(N,))
    rowval = collect(1:N)
    colptr = Vector{Int}(undef,N+1); colptr[1:N] = rowval; colptr[end] = N+1
    nzval  = ones(N)
    return Operator(SparseMatrixCSC(N,N,colptr,rowval,nzval),dims,true)
end
qeye(dims::Dims) = qeye(prod(dims),dims)

"""
    destroy(N)

Generate a quantum harmonic oscillator lowering (annihilation) operator \$\\hat{a}\$ in a truncated Hilbert space of size `N`. Returns a sparse matrix.

# Example
```jldoctest
julia> destroy(4)
4×4 Operator{SparseMatrixCSC{Float64,Int64},1} with dimensions 4
 0.0  1.0  0.0      0.0
 0.0  0.0  1.41421  0.0
 0.0  0.0  0.0      1.73205
 0.0  0.0  0.0      0.0
```
"""
function destroy(N::Integer)
    rowval = collect(1:N-1)
    colptr = Vector{Int}(undef,N+1); colptr[1] = 1; colptr[2:end] = 1:N
    nzval  = [sqrt(i) for i in 1:N-1]
    return Operator(SparseMatrixCSC(N,N,colptr,rowval,nzval),(N,),false)
end

"""
    create(N)

Generate a quantum harmonic oscillator raising (creation) operator \$\\hat{a}^†\$ in a truncated Hilbert space of size `N`. Returns a sparse matrix.

# Example
```jldoctest
julia> create(4)
4×4 Operator{SparseMatrixCSC{Float64,Int64},1} with dimensions 4
 0.0  0.0      0.0      0.0
 1.0  0.0      0.0      0.0
 0.0  1.41421  0.0      0.0
 0.0  0.0      1.73205  0.0
```
"""
function create(N::Integer)
    rowval = collect(2:N)
    colptr = Vector{Int}(undef,N+1); colptr[1:N] = 1:N; colptr[end] = N
    nzval  = [sqrt(i) for i in 1:N-1]
    return Operator(SparseMatrixCSC(N,N,colptr,rowval,nzval),(N,),false)
end

"""
    numberop(N)

Generate a number operator \$\\hat{n}\$ in a Hilbert space of size `N`. Returns a sparse matrix.

# Example
```jldoctest
julia> numberop(4)
4×4 Operator{SparseMatrixCSC{Float64,Int64},1} with dimensions 4
 0.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0
 0.0  0.0  2.0  0.0
 0.0  0.0  0.0  3.0
```
"""
function numberop(N::Integer)
    # "nzval" includes a structural 0 for the [1,1] entry
    rowval = collect(1:N)
    colptr = Vector{Int}(undef,N+1); colptr[1:N] = rowval; colptr[end] = N+1
    nzval  = [float(n) for n = 0:N-1]
    return Operator(SparseMatrixCSC(N,N,colptr,rowval,nzval),(N,),true)
end

"""
    displacementop(N, α)

Generate a quantum harmonic oscillator displacement operator \$\\hat{D}(α)\$ in a truncated Hilbert space of size `N`. Returns a dense matrix.

```math
\\hat{D}(α) = \\exp\\left(α\\hat{a}^† - α^*\\hat{a}\\right)
```

# Example
```jldoctest
julia> displacementop(3,0.5im)
3×3 Operator{Array{Complex{Float64},2},1} with dimensions 3
   0.88262+0.0im            0.0+0.439802im  -0.166001+0.0im
       0.0+0.439802im  0.647859+0.0im             0.0+0.621974im
 -0.166001+0.0im            0.0+0.621974im    0.76524+0.0im
```
"""
function displacementop(N::Integer, α::Number)
    a = full(destroy(N))
    return Operator(exp(α.*a' .- α'.*a),(N,),false)
end

"""
    squeezeop(N, z)

Generate a quantum harmonic oscillator squeeze operator \$\\hat{S}(z)\$ in a truncated Hilbert space of size `N`. Returns a dense matrix.

```math
\\hat{S}(z) = \\exp\\left(\\frac{1}{2}\\left(z^*\\hat{a}^2 - z\\hat{a}^{†2}\\right)\\right)
```

# Example
```jldoctest
julia> squeezeop(3,0.5im)
3×3 Operator{Array{Complex{Float64},2},1} with dimensions 3
 0.938148+0.0im       0.0+0.0im       0.0-0.346234im
      0.0+0.0im       1.0+0.0im       0.0+0.0im
      0.0-0.346234im  0.0+0.0im  0.938148+0.0im
```
"""
function squeezeop(N::Integer, z::Number)
    a = full(destroy(N))
    return Operator(exp(0.5 .* (z'.*a^2 .- z.*a'^2)),(N,),false)
end

"""
    projectorop(N,S)

Generate a projector on the subspaces defined by an integer or a vector/range of integers `S`:

```math
P = \\sum_{i∈S} |i⟩⟨i|.
```

# Example
```jldoctest
julia> projectorop(5,[1,3])
5×5 Operator{SparseMatrixCSC{Float64,Int64},1} with dimensions 5
 0.0  0.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  1.0  0.0
 0.0  0.0  0.0  0.0  0.0
```
"""
function projectorop(N::Integer,S::AbstractVector{<:Integer})
    maximum(S)<N || throw(ArgumentError("a $N-d space cannot be projected on level $(maximum(S))"))
    I = S.+1
    V = ones(length(S))
    return Operator(sparse(I,I,V,N,N),(N,),true)
end
projectorop(N::Integer,S::Integer) = projectorop(N,S:S)

"""
    sylvesterop(N,k,l)

Generate the \$(k,j)^{\textrm{th}}\$ Sylvester generalized Pauli matrix in N-d.
https://en.wikipedia.org/wiki/Generalizations_of_Pauli_matrices
"""
function sylvesterop(N::Integer,k::Integer,j::Integer)
    ωʲ = Complex(cospi(2j/N),sinpi(2j/N))
    rowval = [mod1(i+k,N) for i = 1:N]
    colptr = Vector{Int}(undef,N+1); colptr[1:N] = 1:N; colptr[end] = N+1
    nzval  = [ωʲ^m for m = 0:N-1]
    return Operator(SparseMatrixCSC(N,N,colptr,rowval,nzval),(N,),false)
end

"""
    Sigma1(N)

Generate the N-d shift matrix Σ₁.
"""
function Sigma1(N)
    rowval = circshift(1:N,-1)
    colptr = Vector{Int}(undef,N+1); colptr[1:N] = 1:N; colptr[end] = N+1
    nzval  = ones(N)
    return Operator(SparseMatrixCSC(N,N,colptr,rowval,nzval),(N,),false)
end

"""
    Sigma3(N)

Generate the N-d clock matrix Σ₃.
"""
function Sigma3(N)
    ω = Complex(cospi(2/N),sinpi(2/N))
    rowval = collect(1:N)
    colptr = Vector{undef,Int}(N+1); colptr[1:N] = 1:N; colptr[end] = N+1
    nzval  = [ω^m for m = 0:N-1]
    return Operator(SparseMatrixCSC(N,N,colptr,rowval,nzval),(N,),false)
end
