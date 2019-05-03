using Base: @_inline_meta

function findstate(V::Vector{Ket{T,D}},s::Union{Tuple,Vector}) where {T,D}
    if length(s) != D || any(s .>= dims(V[1]))
        throw(ArgumentError("state $s incompatible with system dimensions $(dims(V[1]))"))
    end
    ti = tensored_sub2ind(dims(V[1]),s)
    j, jval = 1, abs2(V[1][ti]) # initial values
    for i = 2:length(V)
        ival = abs2(V[i][ti])
        if ival > jval; j, jval = i, ival; end
    end
    return j
end

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
