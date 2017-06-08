immutable Result{T<:QuObject,A}
    times::Vector{Float64}
    states::Vector{T}
    evals::Matrix{Complex128}
    #probs::Vector{Matrix{Float64}}
    solver::A
end

function lsolve(L::Liouvillian,ψ₀::Ket,tspan,e_ops,alg;kwargs...)
    dimsmatch(L,ψ₀)
    prob = ODEProblem(L,complex(full(ψ₀)),tspan)
    sol  = solve(prob,alg;dense=false,abstol=1E-8,reltol=1E-6,kwargs...)
    states = Ket.(sol.u,[dims(ψ₀)])
    evals  = calc_expvals(e_ops,states)
    #probs  = calc_probs(states)
    return Result(sol.t,states,evals,sol.alg)
end

function lsolve(L::Liouvillian,ρ₀::Operator,tspan,e_ops,alg;kwargs...)
    dimsmatch(L,ρ₀)
    prob = ODEProblem(L,vec(complex(full(ρ₀))),tspan)
    sol  = solve(prob,alg;dense=false,abstol=1E-8,reltol=1E-6,kwargs...)
    states = Operator.(oper.(sol.u),[dims(ρ₀)])
    evals  = calc_expvals(e_ops,states)
    return Result(sol.t,states,evals,sol.alg)
end

function sesolve(H,ψ₀::Ket,tspan,e_ops=[],alg=Vern8();kwargs...)
    L = SchrodingerEvo(H)
    return lsolve(L,ψ₀,tspan,e_ops,alg;kwargs...)
end

function mesolve(H,C,ρ₀::Operator,tspan,e_ops=[],alg=Tsit5();kwargs...)
    L = LindbladEvo(H,C)
    return lsolve(L,ρ₀,tspan,e_ops,alg;kwargs...)
end

function mesolve(H,C,ψ₀::Ket,tspan,e_ops=[],alg=Tsit5();kwargs...)
    return mesolve(H,C,Operator(ψ₀),tspan,e_ops,alg;kwargs...)
end

function psolve(U::Propagator,ψ₀::Ket,steps,e_ops)
    dimsmatch(U,ψ₀)
    states = Vector{Ket{Vector{Complex{Float64}},length(dims(ψ₀))}}(steps+1)
    states[1] = complex(dense(ψ₀))
    for s = 2:steps+1
        states[s] = U(0.0,states[s-1])
    end
    evals  = calc_expvals(e_ops,states)
    #probs  = calc_probs(states)
    t = collect(linspace(0,steps*U.Δt,steps+1))
    return Result(t,states,evals,:SchrodingerPropSolver)
end

function psolve(U::Propagator,ρ₀::Operator,steps,e_ops)
    dimsmatch(U,ρ₀)
    states = Vector{Operator{Matrix{Complex{Float64}},length(dims(ρ₀))}}(steps+1)
    states[1] = complex(dense(ρ₀))
    for s = 2:steps+1
        states[s] = U(0.0,states[s-1])
    end
    evals  = calc_expvals(e_ops,states)
    #probs  = calc_probs(states)
    t = collect(linspace(0,steps*U.Δt,steps+1))
    return Result(t,states,evals,:SchrodingerPropSolver)
end

function psteady(U::Propagator,ρ₀::Operator,steps,e_ops)
    dimsmatch(U,ρ₀)
    states = Vector{Operator{Matrix{Complex{Float64}},length(dims(ρ₀))}}(2)
    states[1] = complex(dense(ρ₀))
    states[2] = U(0.0,states[1],steps)
    evals  = calc_expvals(e_ops,states)
    #probs  = calc_probs(states)
    return Result([0.0,steps*U.Δt],states,evals,:SchrodingerPropSteady)
end

function calc_expvals(e_ops,states)
    isempty(e_ops) && return Matrix{Complex128}(0,0)
    M = length(e_ops); N = length(states)
    for i=1:M; dimsmatch(states[1],e_ops[i]); end
    expvals = Matrix{Complex128}(N,M)
    for (j,σ) in enumerate(e_ops)
        eval = view(expvals,:,j)
        _expect!(eval, σ, states)
    end
    return expvals
end


tuple_sans_m{n}(m,::Type{Val{n}}) = sorted_setdiff(ntuple(identity,Val{n}),(m,))

calc_probs{T}(ψ::Ket{T,1}) = abs2.(data(ψ))
calc_probs{T,N}(ψ::Ket{T,N},s::Int) = diag(ptrace(ψ,tuple_sans_m(s,Val{N})))
calc_probs{T,N}(ψ::Ket{T,N},out::Vector{Int}) = diag(ptrace(ψ,out))

function calc_probs{T,M}(states::Vector{Ket{T,M}})
    N = length(states)
    S = length(dims(states[1]))
    probs = map(1:S) do s
        d = dims(states[1])[s]
        P = Matrix{Float64}(d,N)
        out = [1:s-1;s+1:N]
        for n = 1:N
            P[:,n] = calc_probs(states[n],s)
        end
        return P.'
    end
    return probs
end

function calc_probs{T}(states::Vector{Ket{T,1}})
    N = length(states)
    d = dims(states[1])[1]
    P = Matrix{Float64}(d,N)
    for n = 1:N
        P[:,n] = calc_probs(states[n])
    end
    return P.'
end

function oper(vecA::AbstractVector,N=isqrt(length(vecA)))
    return reshape(vecA,(N,N))
end

function super(A::AbstractMatrix,B::AbstractMatrix=qeye(size(A,1)).data)
    return B.' ⊗ A
end

# Faster expectation value mathods when evaluating many (see calc_expvals)
function _expect!{T,M}(res,σ::Operator,states::Vector{Ket{T,M}})
    tmp = similar(data(states[1]))
    for (i,ψ) in enumerate(states)
        spmdv_mul!(tmp,data(σ),data(ψ))
        res[i] = dot(data(ψ),tmp)
    end
end

function _expect!{T,M}(res,σ::Operator,states::Vector{Operator{T,M}})
    superσ = super(data(σ))
    tmp = similar(vec(data(states[1])))
    N = length(tmp)
    sqrtNp1 = isqrt(N)+1
    fill!(res, zero(eltype(tmp)))
    for (i,ρ) in enumerate(states)
        spmdv_mul!(tmp,superσ,vec(data(ρ)))
        for n = 1:sqrtNp1:N
            res[i] += tmp[n]
        end
    end
end
