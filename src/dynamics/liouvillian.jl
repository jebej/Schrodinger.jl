immutable Liouvillian{N,F,D} <: AbstractParameterizedFunction{true}
    L₀::SparseMatrixCSC{Complex128,Int}
    Lₙ::NTuple{N,SparseMatrixCSC{Complex128,Int}}
    fₙ::F
    pₙ::NTuple{N,Vector{Float64}}
    dims::SDims{D}
    Liouvillian{N,D}(dims,L₀,Lₙ=(),fₙ=(),pₙ=()) where {N,D} = new{N,typeof(fₙ),D}(L₀,Lₙ,fₙ,pₙ,dims)
end

dimsmatch(L::Liouvillian,A::QuObject) = L.dims==dims(A) || throw(DimensionMismatch("subspace dimensions must match"))

# Liouvillian
function (L::Liouvillian)(t,ψ,dψ)
    A_mul_B!(1.0, L.L₀, ψ, 0.0, dψ)
    applyfun!(dψ, L.Lₙ, L.fₙ, L.pₙ, t, ψ)
end
@inline function applyfun!(dψ,Lₙ,fₙ,pₙ,t,ψ)
    A_mul_B!(first(fₙ)(t,first(pₙ)), first(Lₙ), ψ, 1.0, dψ)
    applyfun!(dψ, tail(Lₙ), tail(fₙ), tail(pₙ), t, ψ)
end
@inline applyfun!(dψ,Lₙ::Tuple{},fₙ,pₙ,t,ψ) = nothing

# Jacobian
function (L::Liouvillian)(::Type{Val{:jac}},t,ψ,J)
    copy!(J, L.L₀)
    applyjac!(J, L.Lₙ, L.fₙ, L.pₙ, t)
end
@inline function applyjac!(J,Lₙ,fₙ,pₙ,t)
    J .+= first(fₙ)(t,first(pₙ)) .* first(Lₙ)
    applyjac!(J, tail(Lₙ), tail(fₙ), tail(pₙ), t)
end
@inline applyjac!(J,Lₙ::Tuple{},fₙ,pₙ,t) = nothing

# Hessian
(p::Liouvillian)(::Type{Val{:hes}},t,ψ,H) = fill!(H, zero(eltype(H)))
