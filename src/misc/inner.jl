const ManyMatrices = AbstractVecOrTuple{T} where T<:AbstractMatrix

@inline inner(A,B) = dot(A,B)

@inline inner(A::NTuple{2,Operator}) = @inbounds inner(A[1],A[2])

@inline inner(A::Tuple{Vararg{Operator}},B::Operator) = inner_cs(data.(A),data(B))

function inner_cs(A,B::AbstractMatrix,s::BitSet=BitSet())
    # calculate |trace(A'*B)|cs (coherent subspaces)
    # this inner product cares only about relative phase between the subspaces contained in s
    x,y = _inner_cs_1(A,B,s)
    return _inner_cs_2(x,y)
end

function inner_cs_grad(Pj,JXj::AbstractMatrix,x,y,s::BitSet=BitSet())
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

function _inner_cs_1(A::AbstractMatrix{T},B::AbstractMatrix{S},s::BitSet) where {T,S}
    m, n = size(A)
    size(B) == (m,n) || throw(DimensionMismatch("matrices must have the same dimensions"))
    x = zero(promote_type(T,S))
    y = Vector{typeof(x)}(undef,n-length(s))
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

function _inner_cs_1(A::ManyMatrices,B::AbstractMatrix,::BitSet)
    # calculate |trace(A'*B)|cs (coherent subspaces)
    # this inner product cares only about relative phase between the subspaces contained in s
    x = inner.(A,(B,))
    return first(x), tail(x)
end

function _inner_cs_1(A::AbstractMatrix,B::ManyMatrices,::BitSet)
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
