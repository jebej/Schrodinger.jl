struct Result{T<:QuObject,A}
    times::Vector{Float64}
    states::Vector{T}
    evals::Matrix{ComplexF64}
    #probs::Vector{Matrix{Float64}}
    solver::A
end

function lsolve(L::Liouvillian,ψ₀::Ket,tspan,e_ops,alg;kwargs...)
    dimsmatch(L,ψ₀)
    prob = ODEProblem(L,complex(full(ψ₀)),tspan)
    sol  = solve(prob,alg;dense=false,abstol=1E-10,reltol=1E-8,kwargs...)
    states = Ket.(sol.u,[dims(ψ₀)])
    evals  = calc_expvals(e_ops,states)
    #probs  = levelprobs(states)
    return Result(sol.t,states,evals,sol.alg)
end

function lsolve(L::Liouvillian,ρ₀::Operator,tspan,e_ops,alg;kwargs...)
    dimsmatch(L,ρ₀)
    prob = ODEProblem(L,vec(complex(full(ρ₀))),tspan)
    sol  = solve(prob,alg;dense=false,abstol=1E-10,reltol=1E-8,kwargs...)
    states = Operator.(oper.(sol.u),[dims(ρ₀)])
    evals  = calc_expvals(e_ops,states)
    return Result(sol.t,states,evals,sol.alg)
end

function sesolve(H,ψ₀::Ket,tspan,e_ops=(),alg=Vern8();kwargs...)
    L = SchrodingerEvo(H)
    return lsolve(L,ψ₀,tspan,e_ops,alg;kwargs...)
end

function mesolve(H,C,ρ₀::Operator,tspan,e_ops=(),alg=Tsit5();kwargs...)
    L = LindbladEvo(H,C)
    return lsolve(L,ρ₀,tspan,e_ops,alg;kwargs...)
end

function mesolve(H,C,ψ₀::Ket,tspan,e_ops=(),alg=Tsit5();kwargs...)
    return mesolve(H,C,Operator(ψ₀),tspan,e_ops,alg;kwargs...)
end

function psolve(U::Propagator,ψ₀::Ket,steps,e_ops)
    dimsmatch(U,ψ₀)
    states = Vector{Ket{Vector{Complex{Float64}},length(dims(ψ₀))}}(steps+1)
    states[1] = complex(dense(ψ₀))
    for s = 2:steps+1
        states[s] = U(states[s-1])
    end
    evals  = calc_expvals(e_ops,states)
    #probs  = levelprobs(states)
    t = collect(linspace(0,steps*U.Δt,steps+1))
    return Result(t,states,evals,:SchrodingerPropSolver)
end

function psolve(U::Propagator,ρ₀::Operator,steps,e_ops)
    dimsmatch(U,ρ₀)
    states = Vector{Operator{Matrix{Complex{Float64}},length(dims(ρ₀))}}(steps+1)
    states[1] = complex(dense(ρ₀))
    for s = 2:steps+1
        states[s] = U(states[s-1])
    end
    evals  = calc_expvals(e_ops,states)
    #probs  = levelprobs(states)
    t = collect(linspace(0,steps*U.Δt,steps+1))
    return Result(t,states,evals,:SchrodingerPropSolver)
end

function psteady(U::Propagator,ρ₀::Operator,steps,e_ops)
    dimsmatch(U,ρ₀)
    states = Vector{Operator{Matrix{Complex{Float64}},length(dims(ρ₀))}}(2)
    states[1] = complex(dense(ρ₀))
    states[2] = U(states[1],steps)
    evals  = calc_expvals(e_ops,states)
    #probs  = levelprobs(states)
    return Result([0.0,steps*U.Δt],states,evals,:SchrodingerPropSteady)
end

calc_expvals(o::Operator,states) = calc_expvals((o,),states)
function calc_expvals(e_ops,states)
    isempty(e_ops) && return Matrix{ComplexF64}(0,0)
    M = length(e_ops); N = length(states)
    expvals = Matrix{ComplexF64}(undef,N,M)
    for (j,σ) in enumerate(e_ops)
        eval = view(expvals,:,j)
        expect!(eval, σ, states)
    end
    return expvals
end
