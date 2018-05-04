import Base: vec, promote_rule, randn

vec(x::RowVector) = x.vec

promote_rule(::Type{SparseMatrixCSC{T1,S1}},::Type{SparseMatrixCSC{T2,S2}}) where {T1,T2,S1,S2} = SparseMatrixCSC{promote_type(T1,T2),promote_type(S1,S2)}
promote_rule(::Type{SparseMatrixCSC{T1,S1}},::Type{<:AbstractMatrix{T2}}) where {T1,S1,T2} = Matrix{promote_type(T1,T2)}
promote_rule(::Type{<:AbstractMatrix{T1}},::Type{SparseMatrixCSC{T2,S2}}) where {T1,T2,S2} = Matrix{promote_type(T1,T2)}

Base.@irrational SQRT_HALF 0.70710678118654752  sqrt(big(0.5))

randn(rng::AbstractRNG,::Type{Complex{T}}) where {T<:AbstractFloat} = Complex{T}(SQRT_HALF*randn(rng,T), SQRT_HALF*randn(rng,T))
