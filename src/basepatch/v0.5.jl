scale!(A::LinAlg.HermOrSym,b::Number) = (scale!(A.data,b);A)

Base.promote_rule{T1,T2,S1,S2}(::Type{SparseMatrixCSC{T1,S1}},::Type{SparseMatrixCSC{T2,S2}}) = SparseMatrixCSC{promote_type(T1,T2),promote_type(S1,S2)}
