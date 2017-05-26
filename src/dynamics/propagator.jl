immutable Propagator{D} <: AbstractParameterizedFunction{false}
    U::Matrix{Complex128}
    Δt::Float64
    #tmp::Vector{Complex128}
    dims::SchroDims{D}
    function (::Type{Propagator{D}}){D}(dims,U,Δt)
        #tmp = Vector{Complex128}(size(L₀,1))
        return new{D}(U,Δt,dims)
    end
end

dimsmatch(U::Propagator,A::QuObject) = U.dims==dims(A) || throw(DimensionMismatch("subspace dimensions must match"))

function (U::Propagator{D}){D}(t,ψ)
    return Ket(U.U*ψ.data,ψ.dims)
end
