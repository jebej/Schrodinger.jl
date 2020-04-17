import Base: BitSet, convert, reverse, print_array, sum, promote_rule
import LinearAlgebra: normalize, BLAS.dotu

const IntCol = Union{AbstractVector{Int},BitSet,Set{Int},NTuple{N,Int} where N}

convert(::Type{BitSet}, r::IntCol) = BitSet(r)

BitSet(elems::Vararg{Int}) = BitSet(elems)

reverse(n::Number) = n

normalize(z::Complex) = z == 0 ? one(z)/abs(one(z)) : z/abs(z)

normalize(z::Real) = one(z)

trace(A::AbstractMatrix) = tr(A)

sum(A::AbstractArray,i::Integer) = Base._sum(A,i)

parseb2(s::AbstractString) = Base.tryparse_internal(Int,s,firstindex(s),lastindex(s),2,true)

const MatrixWrapper{T,S} = Union{Adjoint{T,S}, Transpose{T,S}}
const AllAbstractSparseMatrix{T} = Union{
    AbstractSparseMatrix{T},
    MatrixWrapper{T,<:AbstractSparseMatrix{T}},
    MatrixWrapper{T,<:MatrixWrapper{T,<:AbstractSparseMatrix{T}}}}

promote_rule(::Type{<:AllAbstractSparseMatrix{T1}},::Type{<:AllAbstractSparseMatrix{T2}}) where {T1,T2} = SparseMatrixCSC{promote_type(T1,T2),Int}

promote_rule(::Type{<:AllAbstractSparseMatrix{T1}},::Type{<:AbstractMatrix{T2}}) where {T1,T2} = Matrix{promote_type(T1,T2)}
promote_rule(::Type{<:AbstractMatrix{T1}},::Type{<:AllAbstractSparseMatrix{T2}}) where {T1,T2} = Matrix{promote_type(T1,T2)}

promote_rule(::Type{<:AllAbstractSparseMatrix{T1}},::Type{Matrix{T2}}) where {T1,T2} = Matrix{promote_type(T1,T2)}
promote_rule(::Type{Matrix{T1}},::Type{<:AllAbstractSparseMatrix{T2}}) where {T1,T2} = Matrix{promote_type(T1,T2)}

promote_rule(::Type{Matrix{T1}},::Type{<:MatrixWrapper{T2,<:Matrix{T2}}}) where {T1,T2} = Matrix{promote_type(T1,T2)}
promote_rule(::Type{<:MatrixWrapper{T1,<:Matrix{T1}}},::Type{Matrix{T2}}) where {T1,T2} = Matrix{promote_type(T1,T2)}

promote_rule(::Type{SparseVector{T1,S1}},::Type{SparseVector{T2,S2}}) where {T1,T2,S1,S2} = SparseVector{promote_type(T1,T2),promote_type(S1,S2)}

dotu(x::AbstractVector{T1},y::AbstractVector{T2}) where {T1<:Real,T2<:Real} = dot(x,y)
