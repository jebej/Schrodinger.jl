# Type aliases
@compat const SchroFloat{T<:AbstractFloat} = Union{T, Complex{T}}
@compat const SchroMatrix{T<:SchroFloat} = Union{SparseMatrixCSC{T,Int}, Diagonal{T}, Symmetric{T}, Hermitian{T}, Matrix{T}}
@compat const SchroVector{T<:SchroFloat} = Union{SparseVector{T,Int}, Vector{T}}
@compat const SchroDims{D} = NTuple{D,Int}


# Basic QuObject types
@compat abstract type QuObject end
@compat abstract type QuVector <: QuObject end
@compat abstract type QuMatrix <: QuObject end


"""
    Ket

Ket vector type. This type has two fields, `data` and `dims`, which store the vector data and the subspace dimensions. A `Ket`, like a [`Density`](@ref) matrix or and [`Operator`](@ref) is parameterized by the number of subspaces it lives in. Two different kets must have the same system dimensions in order to be added together.
"""
immutable Ket{T<:SchroVector,D} <: QuVector
    data::T
    dims::SchroDims{D}
    function (::Type{Ket{T,D}}){T<:SchroVector,D}(data,dims)
        prod(dims)==length(data) || throw(ArgumentError("subspace dimensions $dims are not consistent with a vector of length $(length(data))"))
        return new{T,D}(data, dims)
    end
end
Ket{T<:SchroVector,D}(data::T, dims::SchroDims{D}) = Ket{T,D}(data, dims)


"""
    Bra

Bra vector type. The dual vector to the `Ket`
"""
immutable Bra{T<:SchroVector,D} <: QuVector
    data::T
    dims::SchroDims{D}
    function (::Type{Bra{T,D}}){T<:SchroVector,D}(data,dims)
        prod(dims)==length(data) || throw(ArgumentError("subspace dimensions $dims are not consistent with a vector of length $(length(data))"))
        return new{T,D}(data, dims)
    end
end
Bra{T<:SchroVector,D}(data::T, dims::SchroDims{D}) = Bra{T,D}(data, dims)


"""
    Density

Density matrix type. This type has two fields, `data` and `dims`, which store the vector data and the subspace dimensions. A `Density` matrix, like a [`Ket`](@ref) or an [`Operator`](@ref), is parameterized by the number of subspaces it lives in. Two different density matrices must have the same system dimensions in order to be added together. A `Density` matrix is always Hermitian, this is enforced on construction by "hermitianizing" the given matrix.
"""
immutable Density{T<:SchroMatrix,D} <: QuMatrix
    data::T
    dims::SchroDims{D}
    function (::Type{Density{T,D}}){T<:SchroMatrix,D}(data,dims)
        N = checksquare(data)
        prod(dims)==N || throw(ArgumentError("subspace dimensions $dims are not consistent with a matrix of size $N"))
        isapproxhermitian(data) || throw(ArgumentError("a density matrix must be Hermitian"))
        return new{T,D}(data, dims)
    end
end
Density{T<:SchroMatrix,D}(data::T, dims::SchroDims{D}) = Density{T,D}(data, dims)


"""
    Operator

Operator type. An `Operator`, like a [`Ket`](@ref) or a [`Density`](@ref) matrix, is parameterized by the number of subspaces it lives in. Two different density matrices must have the same system dimensions in order to be added together. An `Operator` may or may not be Hermitian.
"""
immutable Operator{T<:SchroMatrix,D} <: QuMatrix
    data::T
    dims::SchroDims{D}
    herm::Bool
    function (::Type{Operator{T,D}}){T<:SchroMatrix,D}(data,dims)
        N = checksquare(data)
        prod(dims)==N || throw(ArgumentError("subspace dimensions $dims are not consistent with a matrix of size $N"))
        return new{T,D}(data, dims, isapproxhermitian(data))
    end
end
Operator{T<:SchroMatrix,D}(data::T, dims::SchroDims{D}) = Operator{T,D}(data, dims)

# A quantum state can be a ket/bra vector or a density matrix
@compat const QuState = Union{Density, Ket, Bra}
