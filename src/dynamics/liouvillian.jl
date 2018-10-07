struct Liouvillian{N,F,D} <: AbstractParameterizedFunction{true}
    L₀::SparseMatrixCSC{ComplexF64,Int}
    Lₙ::NTuple{N,SparseMatrixCSC{ComplexF64,Int}}
    fₙ::F
    pₙ::NTuple{N,Vector{Float64}}
    dims::Dims{D}
    mass_matrix::UniformScaling{Bool}
    analytic::Nothing
end

Liouvillian(dims::Dims{D},L₀,Lₙ=(),fₙ=(),pₙ=()) where {D} =
    Liouvillian{length(Lₙ),typeof(fₙ),D}(L₀,Lₙ,fₙ,pₙ,dims,UniformScaling{Bool}(true),nothing)

dims(L::Liouvillian) = L.dims

dimsmatch(L::Liouvillian,A::QuObject) = dims(L)==dims(A) || throw(DimensionMismatch("subspace dimensions must match"))

# Liouvillian
function (L::Liouvillian)(dψ,ψ,p,t)
    mul!(dψ, L.L₀, ψ, 1, 0)
    applyfun!(dψ, L.Lₙ, L.fₙ, L.pₙ, t, ψ)
end
@inline function applyfun!(dψ,Lₙ,fₙ,pₙ,t,ψ)
    mul!(dψ, first(Lₙ), ψ, first(fₙ)(t,first(pₙ)), 1)
    applyfun!(dψ, tail(Lₙ), tail(fₙ), tail(pₙ), t, ψ)
end
@inline applyfun!(dψ,Lₙ::Tuple{},fₙ,pₙ,t,ψ) = nothing

# Jacobian
function (L::Liouvillian)(::Type{Val{:jac}},J,ψ,p,t)
    copy!(J, L.L₀)
    applyjac!(J, L.Lₙ, L.fₙ, L.pₙ, t)
end
@inline function applyjac!(J,Lₙ,fₙ,pₙ,t)
    J .+= first(fₙ)(t,first(pₙ)) .* first(Lₙ)
    applyjac!(J, tail(Lₙ), tail(fₙ), tail(pₙ), t)
end
@inline applyjac!(J,Lₙ::Tuple{},fₙ,pₙ,t) = nothing
