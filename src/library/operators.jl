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
function qzero(::Type{T}, N::Integer, dims::Dims=(N,)) where {T<:Number}
    return Operator(SparseMatrixCSC{T,Int}(N,N,ones(Int,N+1),Int[],T[]),dims)
end
qzero(::Type{T}, dims::Dims) where {T<:Number} = qzero(T,prod(dims),dims)
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
function qeye(::Type{T}, N::Integer, dims::Dims=(N,)) where {T<:Number}
    return Operator(SparseMatrixCSC{T,Int}(N,N,[1:N+1;],[1:N;],ones(T,N)),dims)
end
qeye(::Type{T}, dims::Dims) where {T<:Number} = qeye(T,prod(dims),dims)
qeye(N::Integer, dims::Dims=(N,)) = qeye(Float64,N,dims)
qeye(dims::Dims) = qeye(Float64,prod(dims),dims)

"""
    destroy(N)

Generate a quantum harmonic oscillator lowering (annihilation) operator ``\\hat{a}`` in a truncated Hilbert space of size `N`. Returns a sparse matrix.

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
function destroy(::Type{T}, N::Integer) where {T<:Number}
    colptr = fill(1,N+1); @inbounds for i in 1:N; colptr[i+1]=i; end
    nzval = T[√(T(i)) for i in 1:N-1]
    return Operator(SparseMatrixCSC{T,Int}(N,N,colptr,[1:N-1;],nzval),(N,))
end
destroy(N::Integer) = destroy(Float64,N)

"""
    create(N)

Generate a quantum harmonic oscillator raising (creation) operator ``\\hat{a}^†`` in a truncated Hilbert space of size `N`. Returns a sparse matrix.

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
function create(::Type{T}, N::Integer) where {T<:Number}
    colptr = fill(N,N+1); @inbounds for i in 1:N; colptr[i]=i; end
    nzval = T[√(T(i)) for i in 1:N-1]
    return Operator(SparseMatrixCSC{T,Int}(N,N,colptr,[2:N;],nzval),(N,))
end
create(N::Integer) = create(Float64,N)

"""
    numberop(N)

Generate a number operator ``\\hat{n}`` in a Hilbert space of size `N`. Returns a sparse matrix.

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
function numberop(::Type{T}, N::Integer) where {T<:Number}
    # nzval includes a structural 0 for the [1,1] entry
    nzval = T[0:N;]
    return Operator(SparseMatrixCSC{T,Int}(N,N,[1:N+1;],[1:N;],nzval),(N,))
end
numberop(N::Integer) = numberop(Float64,N)

"""
    displacementop(N, α)

Generate a quantum harmonic oscillator displacement operator ``\\hat{D}(α)`` in a truncated Hilbert space of size `N`. Returns a dense matrix.

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
    a = Array(destroy(N))
    return Operator(LinearAlgebra.exp!(α.*a' .- α'.*a),(N,))
end

"""
    squeezeop(N, z)

Generate a quantum harmonic oscillator squeeze operator ``\\hat{S}(z)`` in a truncated Hilbert space of size `N`. Returns a dense matrix.

```math
\\hat{S}(z) = \\exp\\left(\\frac{1}{2}\\left(z^*\\hat{a}^2 - z\\hat{a}^{†2}\\right)\\right)
```

# Example
```jldoctest
julia> squeezeop(3,0.5im)
3×3 Operator{Array{Complex{Float64},2},1} with dimensions 3
 0.938148-0.0im       0.0-0.0im       0.0-0.346234im
      0.0+0.0im       1.0+0.0im       0.0+0.0im
      0.0-0.346234im  0.0-0.0im  0.938148-0.0im
```
"""
function squeezeop(N::Integer, z::Number)
    a = Array(destroy(N))
    return Operator(LinearAlgebra.exp!((z'.*a^2 .- z.*a'^2)./2),(N,))
end

"""
    projectorop(N,S)

Generate a projector on the subspaces defined by an integer or a vector/range of integers `S`:

```math
P = ∑_{i∈S} |i⟩⟨i|.
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
    return Operator(sparse(I,I,V,N,N),(N,))
end
projectorop(N::Integer,S::Integer) = projectorop(N,S:S)

"""
    sylvesterop(N,k,l)

Generate the ``(k,j)^\\text{th}`` Sylvester generalized Pauli matrix in N-d.
https://en.wikipedia.org/wiki/Generalizations_of_Pauli_matrices
"""
function sylvesterop(N::Integer,k::Integer,j::Integer)
    ωʲ = Complex(cospi(2j/N),sinpi(2j/N))
    rowval = [mod1(i+k,N) for i in 1:N]
    nzval  = [ωʲ^m for m = 0:N-1]
    return Operator(SparseMatrixCSC(N,N,[1:N+1;],rowval,nzval),(N,))
end

"""
    Sigma1(N)

Generate the N-d shift matrix ``Σ₁``.
"""
function Sigma1(N)
    rowval = circshift(1:N,-1)
    nzval  = ones(N)
    return Operator(SparseMatrixCSC(N,N,[1:N+1;],rowval,nzval),(N,))
end

"""
    Sigma3(N)

Generate the N-d clock matrix ``Σ₃``.
"""
function Sigma3(N)
    ω = Complex(cospi(2/N),sinpi(2/N))
    rowval = [1:N;]
    nzval  = [ω^m for m = 0:N-1]
    return Operator(SparseMatrixCSC(N,N,[1:N+1;],rowval,nzval),(N,))
end
