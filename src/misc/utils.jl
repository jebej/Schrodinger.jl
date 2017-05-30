import Base: reverse, kron, spzeros


size(t::Tuple) = (length(t),)

function sqrfact(n::Integer)
    if 0 <= n <= 20
        return sqrt(factorial(n))
    else # Stirling's formula to avoid overflow for n > 20
        return sqrt(sqrt(2π*n)*(n/e)^n)
    end
end

reverse(n::Number) = n

function expim(A::AbstractMatrix)
    # First decompose A into U*Λ*Uᴴ
    (Λ,U) = eig(A)
    # U must be a complex matrix
    U = complex(U)
    # Calculate the imaginary exponential of each eigenvalue and multiply
    # each column of U to obtain B = U*exp(iΛ)
    n = length(Λ)
    B = similar(U)
    @inbounds for j = 1:n
        a = cis(Λ[j])
        @simd for i = 1:n
            B[i,j] = a * U[i,j]
        end
    end
    # Finally multiply B by Uᴴ to obtain U*exp(iΛ)*Uᴴ = exp(i*A)
    return B*U'
end

function gaussian(x,w)
    return exp(-4ln2*(x/w)^2)
end

spzeros{T<:Number}(a::Type{T},m::Int,n::Int,o::Int) = reshape(spzeros(a,m,n*o),m,n,o)
spzeros{T<:Number}(a::Type{T},m::Int,n::Int,o::Int,p::Int) = reshape(spzeros(a,m,n*o*p),m,n,o,p)

function tensoreye(M::Integer,A::SparseMatrixCSC)
    MN = M*checksquare(A)
    (I,J,V) = findnz(A)
    nz = nnz(A)
    VV = repmat(V,M)
    II = repmat(I,M)
    JJ = repmat(J,M)
    #for m = 0:M-1
    #    lb = m*nz+1; ub = (m+1)*nz; mM = m*M
    #    II[lb:ub] .+= mM; JJ[lb:ub] .+= mM
    #end
    return sparse(II,JJ,VV,MN,MN)
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

function kron{T}(x::Diagonal{T}, y::Diagonal{T})
    valx = x.diag
    valy = y.diag
    nx = length(valx)
    ny = length(valy)
    nz = nx*nx
    valz = Vector{T}(nz)
    @inbounds for i = 1:nx, j = 1:ny
        valz[(i-1)*ny+j] = valx[i] * valy[j]
    end
    return Diagonal(valz)
end

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

function tindex{N}(inds,sysdims::NTuple{N,Int})
    i = inds[N]
    d = 1
    @inbounds for n = reverse(1:N-1)
        d *= sysdims[n+1]
        i += d * (inds[n]-1)
    end
    return i
end
