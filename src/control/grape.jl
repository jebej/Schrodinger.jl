using Optim: Optimizer, Options

abstract type ObjectiveFunction end
abstract type PenaltyFunction end
const IntCol = Union{AbstractVector{Int},IntSet,Set{Int},NTuple{N,Int} where N}

immutable GrapeResult{T<:Operator,S,P<:Optim.MultivariateOptimizationResults}
    Ut::T
    Ui::T
    ui::Matrix{Float64}
    fi::Float64
    Uf::T
    uf::Matrix{Float64}
    ff::Float64
    t::S
    optim_res::P
end

grape(O::ObjectiveFunction,ui::Array,opt::Options) =
grape(O,ui,ConjugateGradient(),opt)

grape(O::ObjectiveFunction,P::PenaltyFunction,ui::Array,opt::Options) =
grape(O,P,ui,ConjugateGradient(),opt)

function grape(O::ObjectiveFunction,ui::Array,
    method::Optimizer=ConjugateGradient(),opt::Options=Options())
    f(u) = objective(O,u)
    g!(fp,u) = gradient!(O,fp,u)
    return grape(f,g!,ui,method,opt)
end

function grape(O::ObjectiveFunction,P::PenaltyFunction,ui::Array,
    method::Optimizer=ConjugateGradient(),opt::Options=Options())
    f(u) = objective(O,u) + objective(P,u)
    g!(fp,u) = (gradient!(O,fp,u); gradient!(P,fp,u))
    return grape(f,g!,ui,method,opt)
end

function grape(f::Function,g!::Function,ui::Array,
    method::Optimizer=ConjugateGradient(),opt::Options=Options())
    # Build OnceDifferentiable object from objective and gradient functions
    fg!(fp,u) = (g!(fp,u); f(u))
    seed = ones(f.O.u_last)
    od = OnceDifferentiable(f,g!,fg!,1.0,similar(seed),seed,copy(seed),[1],[1])
    # Run optimization
    res = optimize(od,vec(ui),method,opt)
    return make_result(f.O, res, ui[:,:])
end

calc_fprops!(O,u) = calc_fprops!(O.U,O.X,O.D,O.V,u,O.δt,O.Hd,O.Hc,O.H,O.u_last)

function calc_fprops!(U,X,D,V,u,δt,Hd,Hc,H,u_last)
    # If the control amplitudes u did not change, return
    u == u_last && return nothing
    # Otherwise calculate each individual propagator
    n = length(U); m = length(Hc)
    for j = 1:n
        copy!(H,Hd)
        for k = 1:m
            H .+= u[(k-1)*n+j].*Hc[k]
        end
        # Calculate U and store eigenvectors and eigenvalues of -δt*H
        expim!(U[j],Hermitian(scale!(H,-δt)),D[j],V[j],X[1])
    end
    # Calculate forward propagators (cumulative product of U)
    copy!(X[1],U[1])
    for i = 2:n
        A_mul_B!(X[i],U[i],X[i-1])
    end
    # Save control amplitudes
    copy!(u_last,u)
    return nothing
end

calc_bprops!(O) = calc_bprops!(O.P,O.U,O.Ut)

function calc_bprops!(P,U,Ut)
    # Calculate backward propagators
    copy!(P[end],Ut)
    for i = length(P)-1:-1:1
        Ac_mul_B!(P[i],U[i+1],P[i+1])
    end
    return nothing
end

function Jmat!(Jkj,Hck,cisDj,Dj,Vj,δt,A)
    # Jₖⱼ = Vⱼ*((Vⱼ'*Hₖ*Vⱼ).*Λⱼ)*Vⱼ'
    # Λⱼ[l,m] = λl≈λm ? -1im*δt*cis(λl) : -δt*(cis(λl)-cis(λm))/(λl-λm)
    Ac_mul_B!(A,Vj,Hck)
    A_mul_B!(Jkj,A,Vj)
    _Jmathermprod!(Jkj,cisDj,Dj,δt)
    A_mul_B!(A,Vj,Jkj)
    A_mul_Bc!(Jkj,A,Vj)
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

function make_result(O,res,ui)
    # Reshape final control amplitude matrix
    uf = reshape(Optim.minimizer(res),size(ui))
    # Calculate final fidelity (should be super cheap)
    ff = 1 - objective(O, Optim.minimizer(res))
    # Grab final operator
    Uf = Operator(O.X[end],O.dims)
    # Calculate initial fidelity
    fi = 1 - objective(O, vec(ui))
    # Grab initial operator
    Ui = Operator(O.X[end],O.dims)
    # Target
    Ut = Operator(O.Ut,O.dims)
    return GrapeResult(Ut,Ui,ui,fi,Uf,uf,ff,O.δt*size(ui,1),res)
end

function plotgrape(res::GrapeResult)
    if :PyPlot ∉ names(Main)
        error("Make sure PyPlot is loaded!")
    else
        plt = Main.PyPlot
    end
    ui = res.ui
    uf = res.uf
    t = linspace(0,res.t,size(uf,1)+1)
    plt.figure()
    plt.step(t,[uf[1:1,:]; uf])
    plt.gca()[:set_prop_cycle](nothing)
    plt.step(t,[ui[1:1,:]; ui],linestyle="dashed")
    plt.legend(["Control $i" for i in 1:size(uf,2)])
    plt.tight_layout(true)
    plt.grid(true)
end
