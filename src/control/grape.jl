using Optim: Optimizer, Options

abstract type ObjectiveFunction end
const IntCol = Union{AbstractVector{Int},IntSet,Set{Int},NTuple{N,Int} where N}

grape(O::ObjectiveFunction,u_init::Array,opt::Options) =
    grape(O,u_init,ConjugateGradient(),opt)

function grape(O::ObjectiveFunction,u_init::Array,
               method::Optimizer=ConjugateGradient(),opt::Options=Options())
    # Build OnceDifferentiable object from ObjectiveFunction object
    f(u) = objective(O,u)
    g!(fp,u) = gradient!(O,fp,u)
    fg!(fp,u) = (g!(fp,u); f(u))
    seed = ones(O.u_last)
    od = OnceDifferentiable(f,g!,fg!,1.0,similar(seed),seed,copy(seed),[1],[1])
    # Run optimization
    res = optimize(od,vec(u_init.'),method,opt)
    # Reshape final control amplitude matrix
    uf = reshape(Optim.minimizer(res),size(u_init,2),size(u_init,1)).'
    return uf, Operator(O.X[end],O.dims), res
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
            H .+= u[(j-1)*m+k].*Hc[k]
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
