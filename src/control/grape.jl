using Optim: AbstractOptimizer, Options, MultivariateOptimizationResults

abstract type ObjectiveFunction end
abstract type PenaltyFunction end

struct GrapeResult{D,T<:AbstractMatrix,P<:MultivariateOptimizationResults}
    Ut::Operator{T,D}
    Ui::Operator{Matrix{ComplexF64},D}
    ui::Matrix{Float64}
    fi::Float64
    Uf::Operator{Matrix{ComplexF64},D}
    uf::Matrix{Float64}
    ff::Float64
    t::Float64
    optim_res::P
end

grape(O::ObjectiveFunction,ui::Array,opt::Options) = grape(O,ui,ConjugateGradient(),opt)

grape(O::ObjectiveFunction,P::PenaltyFunction,ui::Array,opt::Options) = grape(O,P,ui,ConjugateGradient(),opt)

function grape(O::ObjectiveFunction,ui::Array,method::AbstractOptimizer=ConjugateGradient(),opt::Options=Options())
    f(u) = objective(O,u)
    g!(fp,u) = gradient!(O,fp,u)
    return grape(f,g!,ui,method,opt)
end

function grape(O::ObjectiveFunction,P::PenaltyFunction,ui::Array,method::AbstractOptimizer=ConjugateGradient(),opt::Options=Options())
    f(u) = objective(O,u) + objective(P,u)
    g!(fp,u) = (gradient!(O,fp,u); gradient!(P,fp,u))
    return grape(f,g!,ui,method,opt)
end

function grape(f::Function,g!::Function,ui::Array,method::AbstractOptimizer=ConjugateGradient(),opt::Options=Options())
    # Build OnceDifferentiable object from objective and gradient functions
    od = OnceDifferentiable(f,g!,vec(ui))
    # Calculate initial fidelity and initial propagator
    fi = 1 - objective(f.O,ui)
    Ui = current_propagator(f.O)
    # Run optimization
    optim_res = optimize(od,vec(ui),method,opt)
    # Grab optimized control amplitudes
    uf = reshape(Optim.minimizer(optim_res),size(ui,1),size(ui,2))
    # Calculate final fidelity (cheap since propagators are already calculated)
    ff = 1 - objective(f.O,uf)
    # Grab final evolution propagator
    Uf = current_propagator(f.O)
    # Target evolution operator
    Ut = target_propagator(f.O)
    # Return result
    return GrapeResult(Ut,Ui,ui[:,:],fi,Uf,uf,ff,f.O.δt*size(ui,1),optim_res)
end

calc_fprops!(O,u) = calc_fprops!(O.U,O.X,O.D,O.V,u,O.δt,O.Hd,O.Hc,O.H,O.u_last)

function calc_fprops!(U,X,D,V,u,δt,Hd,Hc,H,u_last)
    # If the control amplitudes u did not change, return
    u == u_last && return nothing
    # Otherwise calculate each individual propagator
    n = length(U); m = length(Hc)
    for j = 1:n # loop over timeslices
        H.data .= -δt.*Hd # copy drift (constant) Hamiltonian
        for k = 1:m # loop over controls
            H.data .+= -δt*u[(k-1)*n+j].*Hc[k] # add in time dependent parts
        end
        # Calculate U and store eigenvectors and eigenvalues for this timeslice
        expim!(U[j],D[j],V[j],H,X[1])
    end
    # Calculate forward propagators (cumulative product of U)
    copyto!(X[1],U[1])
    for i = 2:n
        mul!(X[i],U[i],X[i-1])
    end
    # Save control amplitudes
    copyto!(u_last,u)
    return nothing
end

calc_bprops!(O) = calc_bprops!(O.P,O.U,O.Ut)

function calc_bprops!(P::Vector{<:AbstractMatrix},U,Ut)
    # Calculate backward propagators
    for i = length(P)-1:-1:1
        mul!(P[i],adjoint(U[i+1]),P[i+1])
    end
    return nothing
end

function calc_bprops!(P::Vector{<:NTuple{M,AbstractMatrix}},U,Ut) where M
    # Calculate backward propagators
    for i = length(P)-1:-1:1, m = 1:M
        mul!(P[i][m],adjoint(U[i+1]),P[i+1][m])
    end
    return nothing
end

Jmat!(O,k,j) = Jmat!(O.Jkj,O.Hc[k],O.cisDj,O.D[j],O.V[j],O.δt,O.A)

function Jmat!(Jkj,Hck,cisDj,Dj,Vj,δt,A)
    # Jₖⱼ = Vⱼ*((Vⱼ'*Hₖ*Vⱼ).*Λⱼ)*Vⱼ'
    # Λⱼ[l,m] = λl≈λm ? -1im*δt*cis(λl) : -δt*(cis(λl)-cis(λm))/(λl-λm)
    mul!(A,adjoint(Vj),Hck)
    mul!(Jkj,A,Vj)
    _Jmathermprod!(Jkj,cisDj,Dj,δt)
    mul!(A,Vj,Jkj)
    mul!(Jkj,A,adjoint(Vj))
    return nothing
end

function _Jmathermprod!(J,cisDj,Dj,δt)
    for m = 1:length(Dj), l = 1:length(Dj)
        λl, λm = Dj[l], Dj[m]
        cisλl, cisλm = cisDj[l], cisDj[m]
        J[l,m] *= abs(λl-λm)<1E-10 ? -1im*δt*cisλl : -δt*(cisλl-cisλm)/(λl-λm)
    end
    return nothing
end

function current_propagator(O)
    # Return the last propagator calculated
    return Operator(copy(O.X[end]),O.dims)
end

function target_propagator(O)
    # Return the target propagator
    return Operator(O.Ut isa NTuple ? sum(O.Ut) : copy(O.Ut), O.dims)
end

function step_sample(fun,params,tspan,steps)
    t₁, t₂ = tspan
    Δt = (t₂-t₁)/steps
    u = Vector{typeof(fun(t₁,params))}(steps)
    for i = 1:steps
        tᵢ = t₁ + (i-1)*Δt
        u[i] = (fun(tᵢ,params) + fun(tᵢ+0.5Δt,params) + fun(tᵢ+Δt,params))/3
    end
    return u
end
