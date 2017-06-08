# Type aliases
@compat const SFloat{T<:AbstractFloat} = Union{T, Complex{T}}
@compat const SMatrix{T<:SFloat} = Union{SparseMatrixCSC{T,Int}, Diagonal{T}, Symmetric{T}, Hermitian{T}, Matrix{T}}
@compat const SVector{T<:SFloat} = Union{SparseVector{T,Int}, Vector{T}}
@compat const SDims{D} = NTuple{D,Int}

# Base QuObject abstract types
@compat abstract type QuObject end
@compat abstract type QuVector <: QuObject end
@compat abstract type QuMatrix <: QuObject end

# Hermitian indicator type
abstract HermOrNot
immutable Herm <: HermOrNot end
immutable NonHerm <: HermOrNot end


"""
    Ket(x, dims=(length(x),))

Construct a ket state vector from the vector `x`. A vector of length `N` will by default be assumed to be an element of a single Hilbert space of dimension `N`. If the vector is an element of a tensor product of Hilbert spaces, the dimensions can be defined manually by passing a tuple of subspace dimensions `dims`. In that case, `prod(dims)` must equal `length(x)`. By default, the vector is stored in sparse format.

The `Ket` type has two fields, `data` and `dims`, which store the vector data and the subspace dimensions. A `Ket`, like a [`Bra`](@ref) or an [`Operator`](@ref) is parameterized by the number of subspaces it lives in. Two different kets must have the same system dimensions in order to be added together.

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

The `Bra` type has two fields, `data` and `dims`, which store the vector data and the subspace dimensions. A `Bra`, like a [`Ket`](@ref) or an [`Operator`](@ref) is parameterized by the number of subspaces it lives in. Two different kets must have the same system dimensions in order to be added together.

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

"""
    Operator(B, dims=(size(B,1),))

Construct a linear operator from the matrix `B`. An `N`×`N` matrix will by default be assumed to describe an operator that acts on a single Hilbert space of dimension `N`. If the matrix represents a linear operator on a tensor product of Hilbert spaces, the dimensions can be defined manually by passing a tuple of subspace dimensions `dims`. In that case, `prod(dims)` must equal `size(B,1)`.

The `Operator` type has two fields, `data` and `dims`, which store the matrix data and the subspace dimensions. An `Operator`, like a [`Ket`](@ref) or a [`Bra`](@ref), is parameterized by the number of subspaces it lives in. Two different density matrices must have the same system dimensions in order to be added together. An `Operator` may or may not be Hermitian.

# Example
```jldoctest
julia> σ = Operator([0 -im ; im 0])
2×2 Schrodinger.Operator{Schrodinger.Herm,Array{Complex{Float64},2},1} with space dimensions 2:
 0.0+0.0im  0.0-1.0im
 0.0+1.0im  0.0+0.0im
```
"""
immutable Operator{H<:HermOrNot,T<:SMatrix,D} <: QuMatrix
    data::T
    dims::SDims{D}
    function (::Type{Operator{H,T,D}}){H,T<:SMatrix,D}(B,dims)
        N = checksquare(B)
        prod(dims)==N || throw(ArgumentError("subspace dimensions $dims are not consistent with a matrix of size $N"))
        return new{H,T,D}(B, dims)
    end
end
Operator{T<:SMatrix,D}(B::T, dims::SDims{D}=(size(B,1),)) = Operator{isapproxhermitian(B)?Herm:NonHerm,T,D}(B,dims)
Operator(B::AbstractMatrix, dims::SDims=(size(B,1),)) = Operator(float(B),dims)

# Conversion between different QuObjects
Bra(x::Ket) = Bra(conj(x.data),x.dims)
Ket(x::Bra) = Ket(conj(x.data),x.dims)
Operator(x::Ket) = Operator(x.data*x.data',x.dims)
Operator(x::Bra) = Operator(conj(x.data)*x.data.',x.dims)


# Density type remnants
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
