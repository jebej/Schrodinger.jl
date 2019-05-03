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
