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

# the below function is useful to find states at resonance
function findstate(V::Vector{Ket{T,D}},s1::Union{Tuple,Vector},s2::Union{Tuple,Vector}) where {T,D}
    if length(s1) != D || length(s2) != D || any(s1 .>= dims(V[1])) || any(s2 .>= dims(V[1]))
        throw(ArgumentError("state $s incompatible with system dimensions $(dims(V[1]))"))
    end
    ti1 = tensored_sub2ind(dims(V[1]),s1)
    ti2 = tensored_sub2ind(dims(V[1]),s2)
    j1, j1val = 1, zero(abs2(V[1][ti1])) # initial values
    j2, j2val = 1, j1val
    for i = 1:length(V)
        i1val = abs2(V[i][ti1]); i2val = abs2(V[i][ti2])
        if i1val > j1val # found bigger value for state 1
            if abs2(V[j1][ti2]) > j2val # check if current index is better for state 2
                j2, j2val = j1, abs2(V[j1][ti2])
            end
            j1, j1val = i, i1val # now do the swap
            continue
        end
        if i2val > j2val # found bigger value for state 2
            if abs2(V[j2][ti1]) > j1val # check if current index is better for state 1
                j1, j1val = j2, abs2(V[j2][ti1])
            end
            j2, j2val = i, i2val # now do the swap
            continue
        end
    end
    return j1, j2
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
