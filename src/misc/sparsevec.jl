using Base.SparseArrays: nonzeroinds, _spdot

function A_mul_Bf{Tv,Ti}(x::SparseVector{Tv,Ti},y::SparseVector{Tv,Ti},f::Function)
    nnzx = nnz(x)
    nnzy = nnz(y)
    nnzz = nnzx*nnzy # number of nonzeros in new matrix
    I = Vector{Ti}(nnzz) # the indices of nonzeros
    J = Vector{Ti}(nnzz) # the indices of nonzeros
    V  = Vector{Tv}(nnzz) # the values of nonzeros
    @inbounds for i = 1:nnzx, j = 1:nnzy
        this_ind = (i-1)*nnzy+j
        I[this_ind] = x.nzind[i]
        J[this_ind] = y.nzind[j]
        V[this_ind] = x.nzval[i] * f(y.nzval[j])
    end
    return sparse(I,J,V,x.n,y.n)
end
A_mul_Bt{Tv,Ti}(x::SparseVector{Tv,Ti},y::SparseVector{Tv,Ti}) = A_mul_Bf(x,y,identity)
A_mul_Bc{Tv,Ti}(x::SparseVector{Tv,Ti},y::SparseVector{Tv,Ti}) = A_mul_Bf(x,y,conj)

function dotu{Tx<:Number,Ty<:Number}(x::AbstractSparseVector{Tx}, y::AbstractSparseVector{Ty})
    x === y && return sumabs2(x)
    n = length(x)
    length(y) == n || throw(DimensionMismatch())
    xnzind = nonzeroinds(x)
    ynzind = nonzeroinds(y)
    xnzval = nonzeros(x)
    ynzval = nonzeros(y)
    return _spdot(*,1,length(xnzind),xnzind,xnzval,1,length(ynzind),ynzind,ynzval)
end

function dotu{Tx<:Number,Ty<:Number}(x::StridedVector{Tx}, y::AbstractSparseVector{Ty})
    n = length(x)
    length(y) == n || throw(DimensionMismatch())
    nzind = nonzeroinds(y)
    nzval = nonzeros(y)
    s = zero(Tx) * zero(Ty)
    @inbounds for i = 1:length(nzind)
        s += x[nzind[i]] * nzval[i]
    end
    return s
end

function dotu{Tx<:Number,Ty<:Number}(x::AbstractSparseVector{Tx}, y::StridedVector{Ty})
    n = length(y)
    length(x) == n || throw(DimensionMismatch())
    nzind = nonzeroinds(x)
    nzval = nonzeros(x)
    s = zero(Tx) * zero(Ty)
    @inbounds for i = 1:length(nzind)
        s += nzval[i] * y[nzind[i]]
    end
    return s
end
