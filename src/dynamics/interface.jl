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
using OrdinaryDiffEq.DiffEqBase: __solve # for inferability

function lsolve(L::Liouvillian,ψ₀::Ket,tspan,e_ops,alg;kwargs...)
    dimsmatch(L,ψ₀)
    u0 = issparse(ψ₀) ? complex(Array(ψ₀)) : complex(data(ψ₀))
    prob = ODEProblem{true}(LiouvillianODE(L),u0,tspan)
    sol  = __solve(prob,alg;dense=false,abstol=1E-10,reltol=1E-8,kwargs...)
    states = Ket.(convert.(Array,sol.u),(dims(ψ₀),))
    evals  = calc_expvals(e_ops,states)
    #probs  = levelprobs(states)
    return Result(sol.t,states,evals,alg)
end

function lsolve(L::Liouvillian,ρ₀::Operator,tspan,e_ops,alg;kwargs...)
    dimsmatch(L,ρ₀)
    u0 = issparse(ρ₀) ? vec(complex(Array(ρ₀))) :  vec(complex(data(ρ₀)))
    prob = ODEProblem{true}(LiouvillianODE(L),u0,tspan)
    sol  = __solve(prob,alg;dense=false,abstol=1E-10,reltol=1E-8,kwargs...)
    states = Operator.(convert.(Array,unvec.(sol.u)),(dims(ρ₀),))
    evals  = calc_expvals(e_ops,states)
    return Result(sol.t,states,evals,alg)
end

# Steady-state solver

function lsteady(L::Liouvillian,ψ₀::Ket,e_ops;kwargs...)
    # for now input state is used only for dispatch; it could be used as starting guess
    dimsmatch(L,ψ₀)
    sol = solve_steady(L)
    states = [normalize!(Ket(convert(Array,sol),dims(ψ₀)))]
    evals  = calc_expvals(e_ops,states)
    #probs  = levelprobs(states)
    return Result([Inf],states,evals,:SchrodingerLinSolver)
end

function lsteady(L::Liouvillian,ρ₀::Operator,e_ops;kwargs...)
    dimsmatch(L,ρ₀)
    sol = solve_steady(L)
    states = [normalize!(Operator(convert(Array,unvec(sol)),dims(ρ₀)))]
    evals  = calc_expvals(e_ops,states)
    #probs  = levelprobs(states)
    return Result([Inf],states,evals,:SchrodingerLinSolver)
end

function solve_steady(L::Liouvillian)
    # to solve for steady state, compute Liouvillian matrix at t=∞ and solve L*ψ == 0
    L_mat = similar(L.L₀); jac_liouvillian!(L,L_mat,Inf)
    return homogeneous_linsolve(L_mat)
end

homogeneous_linsolve(A::AbstractSparseMatrix) = vec(eigs(A,nev=1,which=:SM,tol=1E-12)[2]::Matrix{eltype(A)})
homogeneous_linsolve(A::AbstractMatrix) = vec(svd(A, full=true).Vt[end,:]')

# Propagator interface

function psolve(U::Propagator,ψ₀::QuObject,steps::Integer,e_ops=())
    dimsmatch(U,ψ₀)
    state1 = complex(dense(ψ₀))
    states = Vector{typeof(state1)}(undef,steps+1); states[1] = state1
    for s in 2:steps+1; states[s] = U(states[s-1]); end
    evals = calc_expvals(e_ops,states)
    #probs  = levelprobs(states)
    return Result(collect(0:U.Δt:steps*U.Δt),states,evals,:SchrodingerPropSolver)
end

function psteady(U::Propagator,ψ₀::QuObject,steps::Integer,e_ops=())
    # same as above but without saving intermediate states
    dimsmatch(U,ψ₀)
    state1 = complex(dense(ψ₀))
    state2 = U(state1,steps)
    states = [state1, state2]
    evals  = calc_expvals(e_ops,states)
    #probs  = levelprobs(states)
    return Result([0.0,steps*U.Δt],states,evals,:SchrodingerPropSteady)
end
