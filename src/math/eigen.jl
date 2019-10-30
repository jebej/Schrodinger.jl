import Compat.LinearAlgebra: eigvals

eigvals(A::QuMatrix, args...; kwargs...) =
    eigvals(A.data, args...; kwargs...)

function eigen(A::QuMatrix, args...; kwargs...)
    D,V = eigen(A.data, args...; kwargs...)
    return D, [Ket(V[:,i],dims(A)) for i = 1:size(V,2)]
end

function eigs(A::Operator, args...; kwargs...)
    D,V = eigs(A.data, args...; which=:SR, kwargs...)
    return D, [Ket(V[:,i],dims(A)) for i = 1:size(V,2)]
end
