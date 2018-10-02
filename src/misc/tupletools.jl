import Base: @_inline_meta

# index in a tensored system (WARNING: subscripts must be 0-indexed, but returns 1-indexed index)
tensored_sub2ind(dims::Dims{N}, inds::Vector{<:Integer}) where {N} =
    (checkbounds(inds,N); @_inline_meta; _tensored_sub2ind(dims, inds, 1, 1))
tensored_sub2ind(dims::Dims{N}, inds::Dims{N}) where {N} =
    (@_inline_meta; _tensored_sub2ind(dims, inds, 1, 1))
_tensored_sub2ind(dims::Dims{N}, inds, D, ind) where {N} =
    (@_inline_meta; _tensored_sub2ind(front(dims), inds, D*dims[N], ind+inds[N]*D))
_tensored_sub2ind(dims::Tuple{Int}, inds, D, ind) = (ind + inds[1]*D)

# give the tensored subscripts corresponding to the index (0-indexed)
tensored_ind2sub(dims::Dims, ind::Integer) =
    (@_inline_meta; next=ind÷dims[end]; (tensored_ind2sub(front(dims), next)..., ind-dims[end]*next))
tensored_ind2sub(dims::Dims{1}, ind::Integer) =
    (@_inline_meta; (ind,))

@inbounds revtuple(t::NTuple{N,Any}) where {N} = ntuple(i->t[N+1-i],Val(N))
@inbounds revinds(t::NTuple{N,Any},ns::Int) where {N} = ntuple(i->ns+1-t[N+1-i],Val(N))
@inbounds gettuple(t1::NTuple,t2::NTuple{N,Any}) where {N} = ntuple(i->t1[t2[i]],Val(N))

setindex(t::Tuple,val,::Val{i}) where {i} = _setindex((),t,val,Val(i))
@inline _setindex(o::Tuple,t::Tuple,val,::Val{i}) where {i} = _setindex((o...,first(t)),tail(t),val,Val(i))
@inline _setindex(o::NTuple{i,Any},t::Tuple,val,::Val{i}) where {i} = (front(o)...,val,t...)

ntuple_sans_m(m,::Val{n}) where {n} = sorted_setdiff(ntuple(identity,Val(n)),(m,))

# Thanks to @mbauman on Discourse for the following
# https://discourse.julialang.org/t/type-stable-difference-of-tuples/3933/4
@inline function sorted_setdiff(t1::Tuple, t2::Tuple)
    if t1[1] == t2[1]
        sorted_setdiff(tail(t1), tail(t2))
    else
        (t1[1], sorted_setdiff(tail(t1), t2)...)
    end
end
@noinline sorted_setdiff(t1::Tuple{}, t2::Tuple) = throw(ArgumentError("duplicate or missing index $(t2[1])"))
sorted_setdiff(t1::Tuple, ::Tuple{}) = t1
sorted_setdiff(::Tuple{}, ::Tuple{}) = ()
