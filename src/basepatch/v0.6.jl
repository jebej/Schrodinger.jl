import Base: kron, promote_rule, randn, vec, dot, vecdot, IntSet
import Base.LinAlg: trace, ctranspose
export sincos

const ComplexF64 = Complex128
const ComplexF32 = Complex64
const adjoint = ctranspose
const qr! = qrfact!

mul!(C::AbstractVecOrMat,A::AbstractVecOrMat,B::AbstractVecOrMat) = A_mul_B!(C,A,B)
mul!(C::AbstractVecOrMat,A::AbstractSparseArray,B::AbstractVecOrMat,a::Number,b::Number) = A_mul_B!(a,A,B,b,C)
dot(A::AbstractMatrix,B::AbstractMatrix) = vecdot(A,B)

vec(x::RowVector) = x.vec

sincos(x::Number) = (sin(x),cos(x))

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

function kron(x::SparseVector{T1,S1}, y::SparseVector{T2,S2}) where {T1,S1,T2,S2}
    nnzx = nnz(x); nnzy = nnz(y)
    nnzz = nnzx*nnzy # number of nonzeros in new vector
    nzind = Vector{promote_type(S1,S2)}(nnzz) # the indices of nonzeros
    nzval = Vector{typeof(one(T1)*one(T2))}(nnzz) # the values of nonzeros
    @inbounds for i = 1:nnzx, j = 1:nnzy
        this_ind = (i-1)*nnzy+j
        nzind[this_ind] = (x.nzind[i]-1)*y.n + y.nzind[j]
        nzval[this_ind] = x.nzval[i] * y.nzval[j]
    end
    return SparseVector(x.n*y.n, nzind, nzval)
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
