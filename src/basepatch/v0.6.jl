import Base: vec, promote_rule

vec(x::RowVector) = x.vec

promote_rule{T1,T2,S1,S2}(::Type{SparseMatrixCSC{T1,S1}},::Type{SparseMatrixCSC{T2,S2}}) = SparseMatrixCSC{promote_type(T1,T2),promote_type(S1,S2)}
promote_rule{T1,S1,T2,S2<:AbstractMatrix{T2}}(::Type{SparseMatrixCSC{T1,S1}},::Type{S2}) = Matrix{promote_type(T1,T2)}
promote_rule{T1,S1,T2,S2<:AbstractMatrix{T2}}(::Type{S2},::Type{SparseMatrixCSC{T1,S1}}) = Matrix{promote_type(T1,T2)}
