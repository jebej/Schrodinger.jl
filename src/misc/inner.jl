function inner{T,S}(A::AbstractMatrix{T},B::AbstractMatrix{S})
    # calculate trace(A'*B) efficiently
    m, n = size(A)
    size(B) == (m,n) || throw(DimensionMismatch("matrices must have the same dimensions"))
    res = zero(promote_type(T,S))
    @inbounds for j = 1:n, i = 1:m
        res += A[i,j]' * B[i,j]
    end
    return res
end

function inner{T1,T2,S1,S2}(A::SparseMatrixCSC{T1,S1},B::SparseMatrixCSC{T2,S2})
    # calculate trace(A'*B) efficiently
    m, n = size(A)
    size(B) == (m,n) || throw(DimensionMismatch("matrices must have the same dimensions"))
    res = zero(promote_type(T1,T2))
    @inbounds for j = 1:n
        for i1 = A.colptr[j]:A.colptr[j+1]-1
            ra = A.rowval[i1]
            for i2 = B.colptr[j]:B.colptr[j+1]-1
                rb = B.rowval[i2]
                if ra < rb
                    # since the rowval of B is larger than that of A, no need to keep checking for equality, go to the nex rowval of A
                    continue
                elseif ra == rb
                    res += A.nzval[i1]' * B.nzval[i2]
                    # done with this row
                    continue
                end
            end
        end
    end
    return res
end

function inner_cs{T,S}(A::AbstractMatrix{T},B::AbstractMatrix{S},s::IntSet)
    # calculate |trace(A'*B)|cs (coherent subspaces)
    # this inner product cares only about relative phase between the subspaces contained in s
    x,y = _inner_cs_1(A,B,s)
    return _inner_cs_2(x,y)
end

function inner_cs_grad{T,S}(Pj::AbstractMatrix{T},JXj::AbstractMatrix{S},x,y,s::IntSet)
    # calculate the partial derivative of |trace(A'*B)|cs (coherent subspaces)
    w,z = _inner_cs_1(Pj,JXj,s)
    res = real(w*x)
    @inbounds for i = 1:length(z)
        res += real(z[i]*y[i])
        res += real(w*normalize(x))*abs(y[i]) + real(z[i]*normalize(y[i]))*abs(x)
        for j = 1:i-1
            res += real(z[i]*normalize(y[i]))*abs(y[j]) + real(z[j]*normalize(y[j]))*abs(y[i])
        end
    end
    return res/_inner_cs_2(x,y)
end

function _inner_cs_1{T,S}(A::AbstractMatrix{T},B::AbstractMatrix{S},s::IntSet)
    m, n = size(A)
    size(B) == (m,n) || throw(DimensionMismatch("matrices must have the same dimensions"))
    x = zero(promote_type(T,S))
    y = Vector{typeof(x)}(n-length(s))
    t = 1
    @inbounds for j = 1:n
        a = zero(x)
        for i = 1:m
            a += A[i,j]' * B[i,j]
        end
        if j in s
            x += a
        else
            y[t] = a
            t += 1
        end
    end
    return x,y
end

function _inner_cs_2(x::Number,y::Vector{<:Number})
    res = abs2(x)
    @inbounds for i = 1:length(y)
        res += abs2(y[i])
        res += 2*abs(y[i])*abs(x)
        for j = 1:i-1
            res += 2*abs(y[i])*abs(y[j])
        end
    end
    return sqrt(res)
end
