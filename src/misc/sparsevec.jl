function A_mul_Bt{Tv,Ti}(x::SparseVector{Tv,Ti},y::SparseVector{Tv,Ti})
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
        V[this_ind] = x.nzval[i] * y.nzval[j]
    end
    return sparse(I,J,V,x.n,y.n)
end

A_mul_Bc{Tv,Ti}(x::SparseVector{Tv,Ti},y::SparseVector{Tv,Ti}) = A_mul_Bt(x,conj(y))
