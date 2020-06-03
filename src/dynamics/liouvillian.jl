struct Liouvillian{T<:Real,M<:AbstractMatrix,N,F,D} <: QuObject
    L₀::M
    Lₙ::NTuple{N,M}
    fₙ::F
    pₙ::NTuple{N,Vector{T}}
    dims::Dims{D}
end

Liouvillian(dims::Dims{D},L₀::M,Lₙ=(),fₙ=(),pₙ=()) where {D,M<:AbstractMatrix} = Liouvillian{real(eltype(M)),M,length(Lₙ),typeof(fₙ),D}(L₀, Lₙ, fₙ, pₙ, dims)

data(L::Liouvillian) = (L.L₀, L.Lₙ, L.fₙ, L.pₙ)

# Liouvillian ODE Type, compatible with AbstractODEFunction interface
using OrdinaryDiffEq.DiffEqBase: AbstractODEFunction
struct LiouvillianODE{LT<:Liouvillian,TD,TJ,TP} <: AbstractODEFunction{true}
    L::LT
    tgrad::TD
    jac::TJ
    jac_prototype::TP
    sparsity::TP
    mass_matrix::UniformScaling{Bool}
end

function LiouvillianODE(L::TL) where {TL<:Liouvillian}
    tgrad = (dT,ψ,p,t) -> tgrad_liouvillian!(L, dT, ψ, t)
    jac = (J,ψ,p,t) -> jac_liouvillian!(L, J, t)
    jac_prot = similar(L.L₀); jac(jac_prot,nothing,nothing,zero(real(eltype(jac_prot))))
    return LiouvillianODE{TL,typeof(tgrad),typeof(jac),typeof(jac_prot)}(L, tgrad, jac, jac_prot, jac_prot, I)
end

(f::LiouvillianODE)(dψ,ψ,p,t) = mul_liouvillian!(f.L, dψ, ψ, t)

# ODE
mul_liouvillian!(L::Liouvillian,dψ,ψ,t) = mul_liouvillian!(dψ, ψ, t, L.L₀, L.Lₙ, L.fₙ, L.pₙ)

function mul_liouvillian!(dψ::AbstractArray,ψ::AbstractArray,t::Real,L₀::AbstractMatrix,Lₙ::Tuple,fₙ::Tuple,pₙ::Tuple)
    mul!(dψ, L₀, ψ, 1, 0)
    applyfun!(dψ, ψ, t, Lₙ, fₙ, pₙ)
end
@inline function applyfun!(dψ::AbstractArray,ψ::AbstractArray,t::Real,Lₙ::Tuple,fₙ::Tuple,pₙ::Tuple)
    mul!(dψ, first(Lₙ), ψ, first(fₙ)(t,first(pₙ)), 1)
    applyfun!(dψ, ψ, t, tail(Lₙ), tail(fₙ), tail(pₙ))
end
@inline applyfun!(dψ::AbstractArray,ψ::AbstractArray,t::Real,Lₙ::Tuple{},fₙ::Tuple{},pₙ::Tuple{}) = nothing

# Jacobian
jac_liouvillian!(L::Liouvillian,J::AbstractMatrix,t::Real) = jac_liouvillian!(J, t, L.L₀, L.Lₙ, L.fₙ, L.pₙ)

function jac_liouvillian!(J::AbstractMatrix,t::Real,L₀::AbstractMatrix,Lₙ::Tuple,fₙ::Tuple,pₙ::Tuple)
    copy!(J, L₀)
    applyjac!(J, t, Lₙ, fₙ, pₙ)
end
@inline function applyjac!(J::AbstractMatrix,t::Real,Lₙ::Tuple,fₙ::Tuple,pₙ::Tuple)
    J += first(fₙ)(t,first(pₙ)) * first(Lₙ)
    applyjac!(J, t, tail(Lₙ), tail(fₙ), tail(pₙ))
end
@inline applyjac!(J::AbstractMatrix,t::Real,Lₙ::Tuple{},fₙ::Tuple{},pₙ::Tuple{}) = nothing

# Time-derivative
tgrad_liouvillian!(L::Liouvillian,dT,ψ,t) = tgrad_liouvillian!(dT, ψ, t, L.L₀, L.Lₙ, L.fₙ, L.pₙ)

function tgrad_liouvillian!(dT::AbstractArray,ψ::AbstractArray,t::Real,L₀::AbstractMatrix,Lₙ::Tuple,fₙ::Tuple,pₙ::Tuple)
    fill!(dT,0)
    apply_tgrad!(dT, ψ, t, Lₙ, fₙ, pₙ)
end
@inline function apply_tgrad!(dT::AbstractArray,ψ::AbstractArray,t::Real,Lₙ::Tuple,fₙ::Tuple,pₙ::Tuple)
    mul!(dT, first(Lₙ), ψ, deriv(t->first(fₙ)(t,first(pₙ)), t), 1)
    apply_tgrad!(dT, ψ, t, tail(Lₙ), tail(fₙ), tail(pₙ))
end
apply_tgrad!(dT::AbstractArray,ψ::AbstractArray,t::Real,Lₙ::Tuple{},fₙ::Tuple{},pₙ::Tuple{}) = nothing

# ForwardDiff doesn't support complex-valued functions so use finite difference in the meantime
function deriv(f::F, x::T) where {F,T<:Real}
    r = sqrt(eps(T))
    h = max(r*abs(x), r)
    return (f(x+h)-f(x-h))/2h
end
