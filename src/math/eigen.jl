import Base: eigfact, eigvals, eig, eigs

eigfact(A::QuMatrix, args...; kwargs...) = eigfact(A.data, args...; kwargs...)

eigvals(A::QuMatrix, args...; kwargs...) = eigvals(A.data, args...; kwargs...)

function eig(A::QuMatrix, args...; kwargs...)
    F = eigfact(A, args...; kwargs...)
    return F.values, [Ket(F.vectors[:,i],dims(A)) for i = 1:size(F.vectors,2)]
end

eigs(A::QuMatrix, args...; kwargs...) = eigs(A.data, args...; kwargs...)
