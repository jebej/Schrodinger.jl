import Base: kron, promote_rule, randn, vec, vecdot

export sincos

vec(x::RowVector) = x.vec

sincos(x::Number) = (sin(x),cos(x))

promote_rule(::Type{SparseMatrixCSC{T1,S1}},::Type{SparseMatrixCSC{T2,S2}}) where {T1,T2,S1,S2} = SparseMatrixCSC{promote_type(T1,T2),promote_type(S1,S2)}
promote_rule(::Type{SparseMatrixCSC{T1,S1}},::Type{<:AbstractMatrix{T2}}) where {T1,S1,T2} = Matrix{promote_type(T1,T2)}
promote_rule(::Type{<:AbstractMatrix{T1}},::Type{SparseMatrixCSC{T2,S2}}) where {T1,T2,S2} = Matrix{promote_type(T1,T2)}

Base.@irrational SQRT_HALF 0.70710678118654752  sqrt(big(0.5))

randn(rng::AbstractRNG,::Type{Complex{T}}) where {T<:AbstractFloat} = Complex{T}(SQRT_HALF*randn(rng,T), SQRT_HALF*randn(rng,T))

function vecdot(A::SparseMatrixCSC{T1,S1},B::SparseMatrixCSC{T2,S2}) where {T1,T2,S1,S2}
    m, n = size(A)
    size(B) == (m,n) || throw(DimensionMismatch("matrices must have the same dimensions"))
    r = zero(promote_type(T1,T2))
    @inbounds for j = 1:n
        ia = A.colptr[j]; ia_nxt = A.colptr[j+1]
        ib = B.colptr[j]; ib_nxt = B.colptr[j+1]
        if ia < ia_nxt && ib < ib_nxt
            ra = A.rowval[ia]; rb = B.rowval[ib]
            while true
                if ra < rb
                    ia += 1
                    ia < ia_nxt || break
                    ra = A.rowval[ia]
                elseif ra > rb
                    ib += 1
                    ib < ib_nxt || break
                    rb = B.rowval[ib]
                else # ra == rb
                    r += A.nzval[ia]' * B.nzval[ib]
                    ia += 1; ib += 1
                    ia < ia_nxt && ib < ib_nxt || break
                    ra = A.rowval[ia]; rb = B.rowval[ib]
                end
            end
        end
    end
    return r
end

function kron{Tv,Ti}(x::SparseVector{Tv,Ti},y::SparseVector{Tv,Ti})
    nnzx = nnz(x)
    nnzy = nnz(y)
    nnzz = nnzx*nnzy # number of nonzeros in new vector
    nzind = Vector{Ti}(nnzz) # the indices of nonzeros
    nzval = Vector{Tv}(nnzz) # the values of nonzeros
    @inbounds for i = 1:nnzx, j = 1:nnzy
        this_ind = (i-1)*nnzy+j
        nzind[this_ind] = (x.nzind[i]-1)*y.n + y.nzind[j]
        nzval[this_ind] = x.nzval[i] * y.nzval[j]
    end
    return SparseVector(x.n*y.n,nzind,nzval)
end

function kron{Tv1,Ti1,Tv2,Ti2}(x::SparseVector{Tv1,Ti1}, y::SparseVector{Tv2,Ti2})
    Tv_res = promote_type(Tv1, Tv2)
    Ti_res = promote_type(Ti1, Ti2)
    A = convert(SparseVector{Tv_res,Ti_res}, x)
    B = convert(SparseVector{Tv_res,Ti_res}, y)
    return kron(A,B)
end

kron(A::SparseVector, B::VecOrMat) = kron(A, sparse(B))
kron(A::VecOrMat, B::SparseVector) = kron(sparse(A), B)

function kron(A::Diagonal{T1}, B::Diagonal{T2}) where {T1<:Number, T2<:Number}
    valA = A.diag; nA = length(valA)
    valB = B.diag; nB = length(valB)
    valC = Vector{typeof(zero(T1)*zero(T2))}(nA*nB)
    @inbounds for i = 1:nA, j = 1:nB
        valC[(i-1)*nB+j] = valA[i] * valB[j]
    end
    return Diagonal(valC)
end
