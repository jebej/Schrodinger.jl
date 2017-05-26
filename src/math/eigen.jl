import Base: eigfact, eigvals, eig, eigs

eigfact(A::QuMatrix, args...; kwargs...) = eigfact(A.data, args...; kwargs...)

eigvals(A::QuMatrix, args...; kwargs...) = eigvals(A.data, args...; kwargs...)

function eig(A::QuMatrix, args...; kwargs...)
    F = eigfact(A, args...; kwargs...)
    F.values, [Ket(F.vectors[:,i],A.dims) for i in 1:size(F.vectors,2)]
end

eigs(A::QuMatrix, args...; kwargs...) = eigs(A.data, args...; kwargs...)
