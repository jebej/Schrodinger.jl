using Schrodinger, Optim
abstract type ObjectiveFunction end
include("normpsu.jl")
include("coherentsubspaces.jl")

function calc_fprops!(U,X,D,V,u,δt,Hd,Hc,H,u_last)
    # If the control amplitudes u did not change, return
    #u == u_last && return nothing
    # Otherwise calculate each individual propagator
    n = length(U); m = length(Hc)
    for j = 1:n
        copy!(H,Hd)
        for k = 1:m
            H .+= u[(j-1)*m+k].*Hc[k]
        end
        # Calculate U and store eigenvectors and eigenvalues of -δt*H
        Schrodinger.expim!(U[j],Hermitian(scale!(H,-δt)),D[j],V[j],X[1])
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

function grape(O::ObjectiveFunction,u_init)
    # Build OnceDifferentiable object from ObjectiveFunction object
    od = OnceDifferentiable(O,(fp,u)->O(Val{:gradient}(),fp,u),zeros(O.u_last))
    # Optimization options
    #opt = Optim.Options(g_tol = 1E-9)
    # Run optimization
    res = optimize(od,vec(u_init.'),ConjugateGradient())
    # Reshape final control amplitude matrix
    uf = reshape(Optim.minimizer(res),size(u_init,2),size(u_init,1)).'
    return uf, O.X[end], res
end

function opt_pihalfx(n)
    Hd = qzero(2)
    Hc = [π*σx, π*σy]
    t = 1
    #u_init = rand(n,1) .- 0.5
    u_init = zeros(n,2)
    Ut = expm(-1im*π*σx/4)
    # Create objective function type
    O = NormPSU(Ut,Hd,Hc,t,n)
    grape(O,u_init)
end

function opt_hadamard(n)
    Hd = σz
    Hc = [σx]
    t = 10.0
    #u_init = rand(n,1) .- 0.5
    u_init = zeros(n,1)
    Ut = Operator([1/√2 1/√2; 1/√2 -1/√2])
    # Create objective function type
    O = NormPSU(Ut,Hd,Hc,t,n)
    grape(O,u_init)
end

function opt_3lvlNOT(n)
    Δ = 2π*(-400) # anharmonicity
    Hd = Δ*Operator(basis(3,2)) # drift Hamiltonian
    Hc = [create(3)/2+destroy(3)/2, im*create(3)/2-im*destroy(3)/2]
    t = 3
    u_init = [-Δ*Schrodinger.gaussianpulse.(linspace(-t/2,t/2,n),[[t/2,t,0,0,pi]]) linspace(-Δ/2,Δ/2,n)]
    Ut = Operator([0 1 0; 1 0 0; 0 0 1]) # 3lvl NOT gate
    # Create objective function type
    O = CoherentSubspaces(Ut,1:2,Hd,Hc,t,n) # care only about computational subspace
    grape(O,u_init)
end

function opt_2Q_QFT(n)
    Hd = 0.5 * (σx⊗σx + σy⊗σy + σz⊗σz)
    Hc = [0.5*σx⊗σ0, 0.5*σy⊗σ0, 0.5*σ0⊗σx, 0.5*σ0⊗σy]
    t = 6.0
    u_init = zeros(n,4)
    #u_init = rand(n,4) .- 0.5
    Ut = Operator([0.5  0.5    0.5  0.5
                   0.5  0.5im -0.5 -0.5im
                   0.5 -0.5    0.5 -0.5
                   0.5 -0.5im -0.5  0.5im],
                  (2,2))
    # Create objective function type
    O = NormPSU(Ut,Hd,Hc,t,n)
    grape(O,u_init)
end

function plotgrape(tf,uf)
    figure()
    step(linspace(0,tf,size(uf,1)+1),[uf[1:1,:]; uf])
    legend(["Control $i" for i in 1:size(uf,2)])
    tight_layout()
    grid()
end
