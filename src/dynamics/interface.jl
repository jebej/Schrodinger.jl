struct Result{T<:QuObject,A}
    times::Vector{Float64}
    states::Vector{T}
    evals::Matrix{ComplexF64}
    #probs::Vector{Matrix{Float64}}
    solver::A
end

# Simple interface: sesolve & mesolve

sesolve(H,ψ₀::Ket,tspan,e_ops=(),alg=Vern8();kwargs...) = lsolve(SchrodingerEvo(H),ψ₀,tspan,e_ops,alg;kwargs...)

mesolve(H,C,ρ₀::Operator,tspan,e_ops=(),alg=Tsit5();kwargs...) = lsolve(LindbladEvo(H,C),ρ₀,tspan,e_ops,alg;kwargs...)

mesolve(H,C,ψ₀::Ket,tspan,e_ops=(),alg=Tsit5();kwargs...) = mesolve(H,C,Operator(ψ₀),tspan,e_ops,alg;kwargs...)

# The functions below require a Liouvillian

function lsolve(L::Liouvillian,ψ₀::Ket,tspan,e_ops,alg;kwargs...)
    dimsmatch(L,ψ₀)
    u0 = complex(Array(ψ₀))
    prob = ODEProblem(LiouvillianODE(L),u0,tspan)
    sol  = solve(prob,alg;dense=false,abstol=1E-10,reltol=1E-8,kwargs...)
    states = Ket.(sol.u::Vector{typeof(u0)},(dims(ψ₀),))
    evals  = calc_expvals(e_ops,states)
    #probs  = levelprobs(states)
    return Result(sol.t::Vector{Float64},states,evals,alg)
end

function lsolve_steady(L::Liouvillian,ψ₀::Ket,e_ops,alg;kwargs...)
    #alg = DynamicSS(odealg;dense=false,abstol=1E-10,reltol=1E-8,kwargs...))
    dimsmatch(L,ψ₀)
    u0 = complex(Array(ψ₀))
    prob = SteadyStateProblem(LiouvillianODE(L),u0)
    sol  = solve(prob,alg;dense=false,abstol=1E-10,reltol=1E-8,kwargs...)
    states = [Ket(sol.u::Vector{typeof(u0)},dims(ψ₀))]
    evals  = calc_expvals(e_ops,states)
    #probs  = levelprobs(states)
    return Result([Inf],states,evals,alg)
end

function lsolve(L::Liouvillian,ρ₀::Operator,tspan,e_ops,alg;kwargs...)
    dimsmatch(L,ρ₀)
    u0 = vec(complex(Array(ρ₀)))
    prob = ODEProblem(LiouvillianODE(L),u0,tspan)
    sol  = solve(prob,alg;dense=false,abstol=1E-10,reltol=1E-8,kwargs...)
    states = Operator.(unvec.(sol.u::Vector{typeof(u0)}),(dims(ρ₀),))
    evals  = calc_expvals(e_ops,states)
    return Result(sol.t::Vector{Float64},states,evals,alg)
end

# Propagator interface

function psolve(U::Propagator,ψ₀::Ket,steps,e_ops)
    dimsmatch(U,ψ₀)
    states = Vector{Ket{Vector{Complex{Float64}},length(dims(ψ₀))}}(undef,steps+1)
    states[1] = complex(dense(ψ₀))
    for s = 2:steps+1
        states[s] = U(states[s-1])
    end
    evals  = calc_expvals(e_ops,states)
    #probs  = levelprobs(states)
    t = collect(LinRange(0,steps*U.Δt,steps+1))
    return Result(t,states,evals,:SchrodingerPropSolver)
end

function psolve(U::Propagator,ρ₀::Operator,steps,e_ops)
    dimsmatch(U,ρ₀)
    states = Vector{Operator{Matrix{Complex{Float64}},length(dims(ρ₀))}}(undef,steps+1)
    states[1] = complex(dense(ρ₀))
    for s = 2:steps+1
        states[s] = U(states[s-1])
    end
    evals  = calc_expvals(e_ops,states)
    #probs  = levelprobs(states)
    t = collect(LinRange(0,steps*U.Δt,steps+1))
    return Result(t,states,evals,:SchrodingerPropSolver)
end

function psteady(U::Propagator,ρ₀::Operator,steps,e_ops)
    dimsmatch(U,ρ₀)
    states = Vector{Operator{Matrix{Complex{Float64}},length(dims(ρ₀))}}(undef,2)
    states[1] = complex(dense(ρ₀))
    states[2] = U(states[1],steps)
    evals  = calc_expvals(e_ops,states)
    #probs  = levelprobs(states)
    return Result([0.0,steps*U.Δt],states,evals,:SchrodingerPropSteady)
end
