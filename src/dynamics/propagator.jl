immutable Propagator{D} <: AbstractParameterizedFunction{false}
    U::Matrix{Complex128}
    Δt::Float64
    dims::SDims{D}
    function (::Type{Propagator{D}}){D}(dims,U,Δt)
        return new{D}(U,Δt,dims)
    end
end

dimsmatch(U::Propagator,A::QuObject) = U.dims==dims(A) || throw(DimensionMismatch("subspace dimensions must match"))

function (U::Propagator{D}){D}(t, ψ::Ket, n::Int=1)
    return Ket(U.U^n*data(ψ),dims(ψ))
end

function (U::Propagator{D}){D}(t, ρ::Density, n::Int=1)
    vecρ = vec(data(ρ))
    d = prod(dims(ρ))
    return Density(reshape(U.U^n*vecρ,(d,d)),dims(ρ))
end
