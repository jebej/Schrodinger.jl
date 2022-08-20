# Type aliases
const SFloat{T<:AbstractFloat} = Union{T, Complex{T}}
const SMatrix{T<:SFloat} = AbstractMatrix{T}
const SVector{T<:SFloat} = AbstractVector{T}

const AbstractVecOrTuple{T} = Union{AbstractVector{<:T},NTuple{N,T} where N}

# Base QuObject abstract types
abstract type QuObject end
abstract type QuVector <: QuObject end
abstract type QuMatrix <: QuObject end

"""
    Ket(x, dims=(length(x),))

Construct a ket state vector from the vector `x`. A vector of length `N` will by default be assumed to be an element of a single Hilbert space of dimension `N`. If the vector is an element of a tensor product of Hilbert spaces, the dimensions can be defined manually by passing a tuple of subspace dimensions `dims`. In that case, `prod(dims)` must equal `length(x)`. By default, the vector is stored in sparse format.

The `Ket` type has two fields, `data` and `dims`, which store the vector data and the subspace dimensions. A `Ket`, like a [`Bra`](@ref) or an [`Operator`](@ref) is parameterized by the number of subspaces it lives in. Two different kets must have the same system dimensions in order to be added together.

It is possible to normalize the ket vector after construction with the `normalize!` function.

# Example
```jldoctest
julia> ψ = normalize!(Ket([1,1]))
2-d Ket{Vector{Float64}, 1} with dimensions 2
0.71∠0°|0⟩ + 0.71∠0°|1⟩
```
"""
struct Ket{T<:SVector,D} <: QuVector
    data::T
    dims::Dims{D}
    function Ket{T,D}(x::AbstractVector, dims::Dims{D}) where {T<:SVector,D}
        checkdims(x,dims)
        return new{T,D}(convert(T,x),dims)
    end
end
Ket(x::T, dims::Dims{D}=size(x)) where {T<:SVector,D} = Ket{T,D}(x,dims)
Ket(x::AbstractVector,dims::Dims=size(x)) = Ket(float.(x),dims)

"""
    Bra(x, dims=(length(x),))

Bra vector type. The dual vector to the `Ket`.

The `Bra` type has two fields, `data` and `dims`, which store the vector data and the subspace dimensions. A `Bra`, like a [`Ket`](@ref) or an [`Operator`](@ref) is parameterized by the number of subspaces it lives in. Two different kets must have the same system dimensions in order to be added together.

It is possible to normalize the bra vector after construction with the `normalize!` function.
"""
struct Bra{T<:SVector,D} <: QuVector
    data::T
    dims::Dims{D}
    function Bra{T,D}(x::AbstractVector,dims::Dims{D}) where {T<:SVector,D}
        checkdims(x,dims)
        return new{T,D}(convert(T,x), dims)
    end
end
Bra(x::T, dims::Dims{D}=size(x)) where {T<:SVector,D} = Bra{T,D}(x,dims)
Bra(x::AbstractVector,dims::Dims=size(x)) = Bra(float.(x),dims)

"""
    Operator(B, dims=(size(B,1),))

Construct a linear operator from the matrix `B`. An `N`×`N` matrix will by default be assumed to describe an operator that acts on a single Hilbert space of dimension `N`. If the matrix represents a linear operator on a tensor product of Hilbert spaces, the dimensions can be defined manually by passing a tuple of subspace dimensions `dims`. In that case, `prod(dims)` must equal `size(B,1)`.

The `Operator` type has two fields, `data` and `dims`, which store the matrix data and the subspace dimensions. An `Operator`, like a [`Ket`](@ref) or a [`Bra`](@ref), is parameterized by the number of subspaces it lives in. Two different density matrices must have the same system dimensions in order to be added together. An `Operator` may or may not be Hermitian.

# Example
```jldoctest
julia> σ = Operator([0 -im ; im 0])
2×2 Operator{Matrix{ComplexF64}, 1} with dimensions 2
 0.0+0.0im  0.0-1.0im
 0.0+1.0im  0.0+0.0im
```
"""
struct Operator{T<:SMatrix,D} <: QuMatrix
    data::T
    dims::Dims{D}
    function Operator{T,D}(B::AbstractMatrix,dims::Dims{D}) where {T<:SMatrix,D}
        checkdims(B,dims)
        return new{T,D}(convert(T,B), dims)
    end
end
Operator(B::T,dims::Dims{D}=(size(B,1),)) where {T<:SMatrix,D} = Operator{T,D}(B,dims)
Operator(B::AbstractMatrix,dims::Dims=(size(B,1),)) = Operator(float.(B),dims)

# Conversion between different QuObjects
Bra(x::Ket) = Bra(conj(x.data),x.dims)
Ket(x::Bra) = Ket(conj(x.data),x.dims)
Operator(x::Ket) = Operator(x.data*x.data',x.dims)
Operator(x::Bra) = Operator(conj(x.data)*tranpose(x.data),x.dims)

# function to check if the dims given are consistent with the object size
function checkdims(x::AbstractVector,dims::Dims{D}) where {D}
    m = length(x)
    if prod(dims) != m
        throw(DimensionMismatch("subspace dimensions $dims are not consistent with a vector of length $m"))
    end
end
function checkdims(A::AbstractMatrix,dims::Dims{D}) where {D}
    m = checksquare(A)
    if prod(dims) != m
        throw(DimensionMismatch("subspace dimensions $dims are not consistent with a matrix of size $m×$m"))
    end
end

# Density type remnants
struct Density{T<:SMatrix,D} <: QuMatrix
    data::T
    dims::Dims{D}
    function Density{T,D}(A::T,dims::Dims{D}) where {T<:SMatrix,D}
        checkdims(A,dims)
        isapproxhermitian(A) || throw(ArgumentError("a density matrix must be Hermitian"))
        return new{T,D}(A, dims)
    end
end
