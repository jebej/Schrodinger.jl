struct Propagator{T,D} <: QuMatrix
    data::Matrix{Complex{T}}
    Δt::T
    dims::Dims{D}
end

function (U::Propagator)(ψ::Ket)
    dimsmatch(U,ψ)
    return Ket(data(U)*data(ψ),dims(ψ))
end

function (U::Propagator)(ψ::Ket, n::Integer)
    dimsmatch(U,ψ)
    return Ket(data(U)^n*data(ψ),dims(ψ))
end

function (U::Propagator)(ρ::Operator)
    dimsmatch(U,ρ)
    d = prod(dims(ρ))
    return Operator(reshape(data(U)*vec(data(ρ)),(d,d)),dims(ρ))
end

function (U::Propagator)(ρ::Operator, n::Integer)
    dimsmatch(U,ρ)
    d = prod(dims(ρ))
    return Operator(reshape(data(U)^n*vec(data(ρ)),(d,d)),dims(ρ))
end

# Convert a Propagator to an Operator
Operator(U::Propagator) = Operator(data(U),dims(U))

# Convert a Liouvillian to a Propagator
function SchrodingerProp(L::Liouvillian{T},tspan,alg=Vern8();kwargs...) where T
    d = prod(dims(L))
    prob = ODEProblem(LiouvillianODE(L),Matrix{Complex{T}}(I,d,d),tspan)
    sol  = __solve(prob,alg;save_start=false,saveat=[tspan[end]],abstol=1E-8,reltol=1E-6,kwargs...)
    return Propagator(sol.u[],(tspan[2]-tspan[1]),dims(L))
end

function LindbladProp(L::Liouvillian{T},tspan,alg=Tsit5();kwargs...) where T
    d = prod(dims(L))^2
    prob = ODEProblem(LiouvillianODE(L),Matrix{Complex{T}}(I,d,d),tspan)
    sol  = __solve(prob,alg;save_start=false,saveat=[tspan[end]],abstol=1E-7,reltol=1E-5,kwargs...)
    return Propagator(sol.u[],(tspan[2]-tspan[1]),dims(L))
end
