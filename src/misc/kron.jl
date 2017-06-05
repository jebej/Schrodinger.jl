import Base: kron

function kron{T1<:Number,T2<:Number}(x::Diagonal{T1}, y::Diagonal{T2})
    T = promote_type(T1, T2)
    valx = x.diag
    valy = y.diag
    nx = length(valx)
    ny = length(valy)
    valz = Vector{T}(nx*ny)
    @inbounds for i = 1:nx, j = 1:ny
        valz[(i-1)*ny+j] = valx[i] * valy[j]
    end
    return Diagonal(valz)
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

function I_kron_A{T}(d::Int, A::AbstractMatrix{T})
    n, m = size(A)
    B = zeros(T,d*n,d*m)
    @inbounds for k = 1:d
        i = (k-1)*n + 1
        j = (k-1)*m + 1
        B[i:i+n-1, j:j+m-1] .= A
    end
    return B
end

function I_kron_A_mul_B!(C,A,B)
    n = checksquare(A)
    n^2 == checksquare(B) == checksquare(C) || throw(DimensionMismatch("invalid dimensions"))
    for j = 1:n
        jj = (j-1)*n+1 : (j-1)*n+n
        for i = 1:n
            ii = (i-1)*n+1 : (i-1)*n+n
            Cij = view(C,ii,jj)
            Bij = view(B,ii,jj)
            A_mul_B!(Cij,A,Bij)
        end
    end
    return C
end

# Not too useful:
function At_kron_I_mul_B!(C,A,B)
    n = checksquare(A)
    n^2 == checksquare(B) == checksquare(C) || throw(DimensionMismatch("invalid dimensions"))
    fill!(C,0)
    kk = 1:n
    for j = 1:n
        jj = (j-1)*n+1 : (j-1)*n+n
        for i = 1:n
            ii = (i-1)*n+1 : (i-1)*n+n
            Cij = view(C,ii,jj)
            for k = 1:n
                Cij .+= A[k,i] .* view(B,kk+(k-1)*n,jj)
            end
        end
    end
    return C
end
