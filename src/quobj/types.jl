# Type aliases
@compat const SFloat{T<:AbstractFloat} = Union{T, Complex{T}}
@compat const SMatrix{T<:SFloat} = Union{SparseMatrixCSC{T,Int}, Diagonal{T}, Symmetric{T}, Hermitian{T}, Matrix{T}}
@compat const SVector{T<:SFloat} = Union{SparseVector{T,Int}, Vector{T}}
@compat const SDims{D} = NTuple{D,Int}


# Basic QuObject types
@compat abstract type QuObject end
@compat abstract type QuVector <: QuObject end
@compat abstract type QuMatrix <: QuObject end


"""
    Ket(x, dims=(length(x),))

Construct a ket state vector from the vector `x`. A vector of length `N` will by default be assumed to be an element of a single Hilbert space of dimension `N`. If the vector is an element of a tensor product of Hilbert spaces, the dimensions can be defined manually by passing a tuple of subspace dimensions `dims`. In that case, `prod(dims)` must equal `length(x)`. By default, the vector is stored in sparse format.

The `Ket` type has two fields, `data` and `dims`, which store the vector data and the subspace dimensions. A `Ket`, like a [`Density`](@ref) matrix or and [`Operator`](@ref) is parameterized by the number of subspaces it lives in. Two different kets must have the same system dimensions in order to be added together.

It is possible to normalize the ket vector after construction with the `normalize!` function.

# Example
```jldoctest
julia> ψ = normalize!(Ket([1,1]))
2-d Schrodinger.Ket{Array{Float64,1},1} with space dimensions 2:
0.71∠0°|0⟩ + 0.71∠0°|1⟩
```
"""
immutable Ket{T<:SVector,D} <: QuVector
    data::T
    dims::SDims{D}
    function (::Type{Ket{T,D}}){T<:SVector,D}(x, dims)
        prod(dims)==length(x) || throw(ArgumentError("subspace dimensions $dims are not consistent with a vector of length $(length(x))"))
        return new{T,D}(x,dims)
    end
end
Ket{T<:SVector,D}(x::T, dims::SDims{D}=(length(x),)) = Ket{T,D}(x,dims)
Ket(x::AbstractVector, dims::SDims=(length(x),)) = Ket(float(x),dims)

"""
    Bra(x, dims=(length(x),))

Bra vector type. The dual vector to the `Ket`.

The `Bra` type has two fields, `data` and `dims`, which store the vector data and the subspace dimensions. A `Bra`, like a [`Density`](@ref) matrix or and [`Operator`](@ref) is parameterized by the number of subspaces it lives in. Two different kets must have the same system dimensions in order to be added together.

It is possible to normalize the bra vector after construction with the `normalize!` function.
"""
immutable Bra{T<:SVector,D} <: QuVector
    data::T
    dims::SDims{D}
    function (::Type{Bra{T,D}}){T<:SVector,D}(x,dims)
        prod(dims)==length(x) || throw(ArgumentError("subspace dimensions $dims are not consistent with a vector of length $(length(x))"))
        return new{T,D}(x, dims)
    end
end
Bra{T<:SVector,D}(x::T, dims::SDims{D}=(length(x),)) = Bra{T,D}(x,dims)
Bra(x::AbstractVector, dims::SDims=(length(x),)) = Bra(float(x),dims)
Bra(x::Ket) = x'
Ket(x::Bra) = x' # Needs to be defined here, after Bra

"""
    Density(A, dims=(size(A,1),))

Construct a density matrix (a.k.a. density operator) from the Hermitian matrix `A`. An `N`×`N` matrix will by default be assumed to describe a single Hilbert space of dimension `N`. If the matrix represents a tensor product of Hilbert spaces, the dimensions can be defined manually by passing a tuple of subspace dimensions `dims`. In that case, `prod(dims)` must equal `size(A,1)`. By default, the matrix is stored in sparse format.

The `Density` type has two fields, `data` and `dims`, which store the matrix data and the subspace dimensions. A `Density` matrix, like a [`Ket`](@ref) or an [`Operator`](@ref), is parameterized by the number of subspaces it lives in. Two different density matrices must have the same system dimensions in order to be added together. A `Density` matrix is always Hermitian, this is enforced on construction by "hermitianizing" the given matrix.

It is possible to normalize the density matrix after construction with the `normalize!` function.

# Example
```jldoctest
julia> A = [1 5 2; 5 1 0; 2 0 2]
3×3 Array{Int64,2}:
 1  5  2
 5  1  0
 2  0  2
julia> σ = normalize!(Density(A))
3×3 Schrodinger.Density{Array{Float64,2},1} with space dimensions 3:
 0.25  1.25  0.5
 1.25  0.25  0.0
 0.5   0.0   0.5
julia> trace(σ)
1.0
```
"""
immutable Density{T<:SMatrix,D} <: QuMatrix
    data::T
    dims::SDims{D}
    function (::Type{Density{T,D}}){T<:SMatrix,D}(A,dims)
        N = checksquare(A)
        prod(dims)==N || throw(ArgumentError("subspace dimensions $dims are not consistent with a matrix of size $N"))
        isapproxhermitian(A) || throw(ArgumentError("a density matrix must be Hermitian"))
        return new{T,D}(A, dims)
    end
end
Density{T<:SMatrix,D}(A::T, dims::SDims{D}=(size(A,1),)) = Density{T,D}(A,dims)
Density(A::AbstractMatrix, dims::SDims=(size(A,1),)) = Density(float(A),dims)
Density(x::Ket) = Density(x.data*x.data',x.dims)
Density(x::Bra) = Density(conj(x.data)*x.data.',x.dims)


"""
    Operator(B, dims=(size(B,1),))

Construct a linear operator from the matrix `B`. An `N`×`N` matrix will by default be assumed to describe an operator that acts on a single Hilbert space of dimension `N`. If the matrix represents a linear operator on a tensor product of Hilbert spaces, the dimensions can be defined manually by passing a tuple of subspace dimensions `dims`. In that case, `prod(dims)` must equal `size(B,1)`.

The `Operator` type has two fields, `data` and `dims`, which store the matrix data and the subspace dimensions. An `Operator`, like a [`Ket`](@ref) or a [`Density`](@ref) matrix, is parameterized by the number of subspaces it lives in. Two different density matrices must have the same system dimensions in order to be added together. An `Operator` may or may not be Hermitian.

# Example
```jldoctest
julia> σ = Operator([0 -im ; im 0])
2×2 Schrodinger.Operator{Array{Complex{Float64},2},1} with space dimensions 2:
 0.0+0.0im  0.0-1.0im
 0.0+1.0im  0.0+0.0im
```
"""
immutable Operator{T<:SMatrix,D} <: QuMatrix
    data::T
    dims::SDims{D}
    herm::Bool
    function (::Type{Operator{T,D}}){T<:SMatrix,D}(B,dims)
        N = checksquare(B)
        prod(dims)==N || throw(ArgumentError("subspace dimensions $dims are not consistent with a matrix of size $N"))
        return new{T,D}(B, dims, isapproxhermitian(B))
    end
end
Operator{T<:SMatrix,D}(B::T, dims::SDims{D}=(size(B,1),)) = Operator{T,D}(B,dims)
Operator(B::AbstractMatrix, dims::SDims=(size(B,1),)) = Operator(float(B),dims)
Operator(x::Ket) = Operator(x.data*x.data',x.dims)
Operator(x::Bra) = Operator(conj(x.data)*x.data.',x.dims)

# A quantum state can be a ket/bra vector or a density matrix
@compat const QuState = Union{Ket, Bra, Density}
