import Base: promote_rule
import LinearAlgebra: BLAS.dotu

promote_rule(::Type{SparseMatrixCSC{T1,S1}},::Type{SparseMatrixCSC{T2,S2}}) where {T1,T2,S1,S2} = SparseMatrixCSC{promote_type(T1,T2),promote_type(S1,S2)}
promote_rule(::Type{SparseMatrixCSC{T1,S1}},::Type{<:AbstractMatrix{T2}}) where {T1,S1,T2} = Matrix{promote_type(T1,T2)}
promote_rule(::Type{<:AbstractMatrix{T1}},::Type{SparseMatrixCSC{T2,S2}}) where {T1,T2,S2} = Matrix{promote_type(T1,T2)}

promote_rule(::Type{SparseVector{T1,S1}},::Type{SparseVector{T2,S2}}) where {T1,T2,S1,S2} = SparseVector{promote_type(T1,T2),promote_type(S1,S2)}

dotu(x::AbstractVector{T1},y::AbstractVector{T2}) where {T1<:Real,T2<:Real} = dot(x,y)
