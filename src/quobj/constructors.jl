# Convenience contrcutors for QuObjects, these should be used, and not the raw constructor

"""
    ket(x, dims=(length(x),), dense=false)

Construct a ket state vector from the vector `x`. A vector of length `N` will by default be assumed to be an element of a single Hilbert space of dimension `N`. If the vector is an element of a tensor product of Hilbert spaces, the dimensions can be defined manually by passing a tuple of subspace dimensions `dims`. In that case, `prod(dims)` must equal `length(a)`. By default, the vector is stored in sparse format.

It is possible to normalize the ket vector after construction with the `normalize!` function.

# Example
```jldoctest
julia> ψ = normalize!(ket([1,1]))
2-d Schrodinger.Ket{SparseVector{Float64,Int64},1} with space dimensions 2:
0.71∠0°|0⟩ + 0.71∠0°|1⟩
```
"""
function ket(x::SchroVector, dims::SchroDims=(length(x),), dense::Bool=false)
    if dense
        return Ket(full(x), dims)
    else
        return Ket(sparse(x), dims)
    end
end
function ket(x::AbstractVector, dims::SchroDims=(length(x),), dense::Bool=false)
    return ket(float.(x), dims, dense)
end

function bra(x::SchroVector, dims::SchroDims=(length(x),), dense::Bool=false)
    if dense
        return Bra(full(x), dims)
    else
        return Bra(sparse(x), dims)
    end
end
function bra(x::AbstractVector, dims::SchroDims=(length(x),), dense::Bool=false)
    return bra(float.(x), dims, dense)
end


"""
    density(A, dims=(size(A,1),), dense=false)

Construct a density matrix (or density operator) from the Hermitian matrix `A`. An `N`×`N` matrix will by default be assumed to describe a single Hilbert space of dimension `N`. If the matrix represents a tensor product of Hilbert spaces, the dimensions can be defined manually by passing a tuple of subspace dimensions `dims`. In that case, `prod(dims)` must equal `size(A,1)`. By default, the matrix is stored in sparse format.

It is possible to normalize the density matrix after construction with the `normalize!` function.

# Example
```jldoctest
julia> A = [1 5 2; 5 1 0; 2 0 2]
3×3 Array{Int64,2}:
 1  5  2
 5  1  0
 2  0  2
julia> σ = normalize!(density(A))
3×3 Schrodinger.Density{SparseMatrixCSC{Float64,Int64},1} with space dimensions 3:
 0.25  1.25  0.5
 1.25  0.25  0.0
 0.5   0.0   0.5
julia> trace(σ)
1.0
```
"""
function density(A::SchroMatrix, dims::SchroDims=(size(A,1),), dense::Bool=false)
    if dense
        return Density(full(A), dims)
    else
        return Density(sparse(A), dims)
    end
end
function density(A::AbstractMatrix, dims::SchroDims=(size(A,1),), dense::Bool=false)
    return density(float.(A), dims, dense)
end
density(x::Ket, dense::Bool=!issparse(x.data)) = density(x.data*x.data', x.dims, dense)
density(x::Bra, dense::Bool=!issparse(x.data)) = density(conj(x.data)*x.data.', x.dims, dense)


"""
    operator(B, dims=(size(B,1),), dense=false)

Construct a linear operator from the matrix `B`. An `N`×`N` matrix will by default be assumed to describe an operator that acts on a single Hilbert space of dimension `N`. If the matrix represents a linear operator on a tensor product of Hilbert spaces, the dimensions can be defined manually by passing a tuple of subspace dimensions `dims`. In that case, `prod(dims)` must equal `size(B,1)`. By default, the matrix is stored in sparse format.

# Example
```jldoctest
julia> σ = operator([0 -im ; im 0])
2×2 Schrodinger.Operator{SparseMatrixCSC{Complex{Float64},Int64},1} with space dimensions 2:
 0.0+0.0im  0.0-1.0im
 0.0+1.0im  0.0+0.0im
```
"""
function operator(B::SchroMatrix, dims::SchroDims=(size(B,1),), dense::Bool=false)
    if dense
        B = full(B)
        isreal(B) && (B = real.(B))
        isdiag(B)      && return Operator(Diagonal(B), dims)
        issymmetric(B) && return Operator(Symmetric(B), dims)
        ishermitian(B) && return Operator(Hermitian(B), dims)
        return Operator(B, dims)
    else
        return Operator(sparse(B), dims)
    end
end
function operator(B::AbstractMatrix, dims::SchroDims=(size(B,1),), dense::Bool=false)
    return operator(float.(B), dims, dense)
end
operator(x::Ket, dense::Bool=!issparse(x.data)) = operator(x.data*x.data', x.dims, dense)
operator(x::Bra, dense::Bool=!issparse(x.data)) = operator(conj(x.data)*x.data.', x.dims, dense)
