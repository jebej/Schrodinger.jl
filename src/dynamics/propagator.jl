struct Propagator{T,D} <: AbstractParameterizedFunction{false}
    U::Matrix{Complex{T}}
    Δt::T
    dims::Dims{D}
end

dimsmatch(U::Propagator,A::QuObject) = U.dims==dims(A) || throw(DimensionMismatch("subspace dimensions must match"))

function (U::Propagator)(ψ::Ket)
    dimsmatch(U,ψ)
    return Ket(U.U*data(ψ),dims(ψ))
end

function (U::Propagator)(ψ::Ket, n::Integer)
    dimsmatch(U,ψ)
    return Ket(U.U^n*data(ψ),dims(ψ))
end

function (U::Propagator)(ρ::Operator)
    dimsmatch(U,ρ)
    vecρ = vec(data(ρ))
    d = prod(dims(ρ))
    return Operator(reshape(U.U*vecρ,(d,d)),dims(ρ))
end

function (U::Propagator)(ρ::Operator, n::Integer)
    dimsmatch(U,ρ)
    vecρ = vec(data(ρ))
    d = prod(dims(ρ))
    return Operator(reshape(U.U^n*vecρ,(d,d)),dims(ρ))
end

# Convert a Propagator to an Operator
Operator(U::Propagator) = Operator(U.U,U.dims)

# Convert a Liouvillian to a Propagator
function SchrodingerProp(L::Liouvillian,tspan,alg=Vern8();kwargs...)
    prob = ODEProblem(L,eye(ComplexF64,prod(dims(L))),tspan)
    sol  = solve(prob,alg;save_start=false,saveat=tspan[end],abstol=1E-8,reltol=1E-6,kwargs...)
    return Propagator(sol.u[end],float(tspan[2]-tspan[1]),dims(L))
end

function LindbladProp(L::Liouvillian,tspan,alg=Tsit5();kwargs...)
    prob = ODEProblem(L,eye(ComplexF64,prod(dims(L))^2),tspan)
    sol  = solve(prob,alg;save_start=false,saveat=tspan[end],abstol=1E-7,reltol=1E-5,kwargs...)
    return Propagator(sol.u[end],float(tspan[2]-tspan[1]),dims(L))
end
