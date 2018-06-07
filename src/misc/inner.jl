const ManyMatrices = Union{Vector{<:AbstractMatrix},NTuple{M,<:AbstractMatrix} where M}

# Remove the two inner functions below (dense and sparse) for 0.7

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

function inner(A::SparseMatrixCSC{T1,S1},B::SparseMatrixCSC{T2,S2}) where {T1,T2,S1,S2}
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

function inner_cs(A,B::AbstractMatrix,s::IntSet)
    # calculate |trace(A'*B)|cs (coherent subspaces)
    # this inner product cares only about relative phase between the subspaces contained in s
    x,y = _inner_cs_1(A,B,s)
    return _inner_cs_2(x,y)
end

function inner_cs_grad(Pj,JXj::AbstractMatrix,x,y,s::IntSet)
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

# sub-routines

function _inner_cs_1(A::AbstractMatrix{T},B::AbstractMatrix{S},s::IntSet) where {T,S}
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
    return x, y
end

function _inner_cs_1(A::ManyMatrices,B::AbstractMatrix,::IntSet)
    # calculate |trace(A'*B)|cs (coherent subspaces)
    # this inner product cares only about relative phase between the subspaces contained in s
    x = inner.(A,(B,))
    return first(x), tail(x)
end

function _inner_cs_1(A::AbstractMatrix,B::ManyMatrices,::IntSet)
    # calculate |trace(A'*B)|cs (coherent subspaces)
    # this inner product cares only about relative phase between the subspaces contained in s
    x = inner.((A,),B)
    return first(x), tail(x)
end

function _inner_cs_2(x::Number,y::Union{Vector{<:Number},NTuple{M,<:Number} where M})
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

f_cs_2(x::Number,y::Union{Vector{<:Number},NTuple{M,<:Number} where M}) = sum(abs,y)+abs(x)
