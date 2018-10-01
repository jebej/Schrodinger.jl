import Compat.LinearAlgebra: eigfact, eigvals, eig, eigs

eigfact(A::QuMatrix, args...; kwargs...) = eigfact(A.data, args...; kwargs...)

eigvals(A::QuMatrix, args...; kwargs...) = eigvals(A.data, args...; kwargs...)

function eig(A::QuMatrix, args...; kwargs...)
    F = eigfact(A, args...; kwargs...)
    return F.values, [Ket(F.vectors[:,i],dims(A)) for i = 1:size(F.vectors,2)]
end

function eig(A::Operator{<:SparseMatrixCSC,N}, args...; kwargs...) where N
    D,V,_ = eigs(A.data, args...; which=:SR, kwargs...)
    return return D, [Ket(V[:,i],dims(A)) for i = 1:size(V,2)]
end
