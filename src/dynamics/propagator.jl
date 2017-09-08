immutable Propagator{D} <: AbstractParameterizedFunction{false}
    U::Matrix{Complex128}
    Δt::Float64
    dims::SDims{D}
    function (::Type{Propagator{D}}){D}(dims,U,Δt)
        return new{D}(U,Δt,dims)
    end
end

dimsmatch(U::Propagator,A::QuObject) = U.dims==dims(A) || throw(DimensionMismatch("subspace dimensions must match"))

function (U::Propagator{D}){D}(ψ::Ket)
    dimsmatch(U,ψ)
    return Ket(U.U*data(ψ),dims(ψ))
end

function (U::Propagator{D}){D}(ψ::Ket, n::Integer)
    dimsmatch(U,ψ)
    return Ket(U.U^n*data(ψ),dims(ψ))
end

function (U::Propagator{D}){D}(ρ::Operator)
    dimsmatch(U,ρ)
    vecρ = vec(data(ρ))
    d = prod(dims(ρ))
    return Operator(reshape(U.U*vecρ,(d,d)),dims(ρ))
end

function (U::Propagator{D}){D}(ρ::Operator, n::Integer)
    dimsmatch(U,ρ)
    vecρ = vec(data(ρ))
    d = prod(dims(ρ))
    return Operator(reshape(U.U^n*vecρ,(d,d)),dims(ρ))
end

# method that takes in t for DifferentialEquation's steady state solver
(U::Propagator{D}){D}(t::Real, ψ::Ket) = U(ψ)
