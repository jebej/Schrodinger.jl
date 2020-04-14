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

promote_rule(::Type{SparseMatrixCSC{T1,S1}},::Type{SparseMatrixCSC{T2,S2}}) where {T1,T2,S1,S2} = SparseMatrixCSC{promote_type(T1,T2),promote_type(S1,S2)}
promote_rule(::Type{SparseMatrixCSC{T1,S1}},::Type{<:AbstractMatrix{T2}}) where {T1,S1,T2} = Matrix{promote_type(T1,T2)}
promote_rule(::Type{<:AbstractMatrix{T1}},::Type{SparseMatrixCSC{T2,S2}}) where {T1,T2,S2} = Matrix{promote_type(T1,T2)}

promote_rule(::Type{SparseVector{T1,S1}},::Type{SparseVector{T2,S2}}) where {T1,T2,S1,S2} = SparseVector{promote_type(T1,T2),promote_type(S1,S2)}

dotu(x::AbstractVector{T1},y::AbstractVector{T2}) where {T1<:Real,T2<:Real} = dot(x,y)
