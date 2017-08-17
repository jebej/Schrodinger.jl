using Schrodinger, Optim

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
        # Calculate U and store eigenvectors and eigenvalues
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
        if abs(λl-λm) < 1E-10
            J[l,m] *= -1im*δt*cisλl
        else
            J[l,m] *= -δt*(cisλl-cisλm)/(λl-λm)
        end
    end
    return nothing
end

function fidelity!(u,δt,Ut,Hd,Hc,H,U,X,D,V,u_last)
    # Calculate forward propagators
    calc_fprops!(U,X,D,V,u,δt,Hd,Hc,H,u_last)
    Uf = X[end] # Full propagator
    # Calculate fidelity error:
    # fₑ = 1 - Φ where Φ = |⟨Ut,Uf⟩|²/N²
    Φ = abs2(inner(Ut,Uf))/size(Ut,1)^2
    return 1 - Φ
end

function fidelityprime!(fp,u,δt,Ut,Hd,Hc,H,A,U,X,D,V,P,u_last)
    # Calculate forward and backward propagators
    calc_fprops!(U,X,D,V,u,δt,Hd,Hc,H,u_last)
    calc_bprops!(P,U,Ut)
    n = length(U); m = length(Hc)
    # Calculate exact derivative of fidelity error function:
    # ∂Φ/∂uₖⱼ  = ⟨Pⱼ,∂Uⱼ/∂uₖⱼ*Xⱼ₋₁⟩⟨Xⱼ,Pⱼ⟩/N² + c.c.
    # ∂fₑ/∂uₖⱼ = -∂Φ/∂uₖⱼ
    #          = -2*Re(⟨Pⱼ,Jₖⱼ*Xⱼ₋₁⟩⟨Uf,Ut⟩)/N² where Jₖⱼ = ∂Uⱼ/∂uₖⱼ
    Jkj = similar(A)
    cisDj = Vector{Complex128}(D[1])
    a = inner(X[end],Ut)
    for j = 1:n
        cisDj .= cis.(D[j])
        for k = 1:m
            Jmat!(Jkj,Hc[k],cisDj,D[j],V[j],δt,A)
            j==1 ? copy!(A,Jkj) : A_mul_B!(A,Jkj,X[j-1])
            fp[(j-1)*m+k] = -real(inner(P[j],A)*a)
        end
    end
    scale!(fp,2/size(Ut,1)^2)
    return fp
end

function fidelityprimeapprox!(fp,u,δt,Ut,Hd,Hc,H,A,U,X,D,V,P,u_last)
    # Calculate forward and backward propagators
    calc_fprops!(U,X,D,V,u,δt,Hd,Hc,H,u_last)
    calc_bprops!(P,U,Ut)
    n = length(U); m = length(Hc)
    # Calculate approximate derivative of fidelity error function:
    # ∂Φ/∂uₖⱼ  = ⟨Pⱼ,∂Uⱼ/∂uₖⱼ*Xⱼ₋₁⟩⟨Xⱼ,Pⱼ⟩/N² + c.c.
    # ∂fₑ/∂uₖⱼ = -∂Φ/∂uₖⱼ
    #          ≈ 2*Re(⟨Pⱼ,iδt*Hₖ*Xⱼ⟩⟨Uf,Ut⟩)/N²
    a = inner(X[end],Ut)
    for j = 1:n
        for k = 1:m
            A_mul_B!(A,Hc[k],X[j])
            fp[(j-1)*m+k] = real(inner(P[j],scale!(A,1im*δt))*aj)
        end
    end
    scale!(fp,2/size(Ut,1)^2)
    return fp
end

function gen_opt_fun(Ut::Operator,Hd::Operator,Hc::Vector{<:Operator},t::Real,n::Integer)
    N = prod(dims(Hd))
    m = length(Hc)
    # Make sure we pass dense operators
    Ut_d = full(Ut)
    Hd_d = full(Hd)
    Hc_d = full.(Hc)
    # Generate cache for various objects
    H = promote_type(typeof(Hd_d),eltype(Hc_d))(N,N)
    A = Matrix{Complex128}(N,N)
    U = [similar(A) for i=1:n]
    X = deepcopy(U)
    P = deepcopy(U)
    # For the exact derivative we need to store eigenvectors and eigenvalues
    V = [similar(H) for i=1:n]
    D = [Vector{Float64}(N) for i=1:n] # eigenvalues are always real
    # Storage for last control ampitudes, NaN for first run
    u_last = fill(NaN64,m*n)
    # Create optimization function object
    f = (u) -> fidelity!(u,t/n,Ut_d,Hd_d,Hc_d,H,U,X,D,V,u_last)
    g! = (fp,u) -> fidelityprime!(fp,u,t/n,Ut_d,Hd_d,Hc_d,H,A,U,X,D,V,P,u_last)
    #g! = (fp,u) -> fidelityprimeapprox!(fp,u,t/n,Ut_d,Hd_d,Hc_d,H,A,U,X,D,V,P,u_last)
    # Return a OnceDifferentiable object with appropriate seed
    return OnceDifferentiable(f,g!,zeros(m*n)),X[end]
end


function grape(Ut::Operator,Hd::Operator,Hc::Vector{<:Operator},u_init,t::Real,n::Integer)
    m = length(Hc)
    n == size(u_init,1) || throw(ArgumentError("control amplitude matrix not consistent with number of timesteps"))
    m == size(u_init,2) || throw(ArgumentError("control amplitude matrix not consistent with number of control hamiltonians"))
    # Generate OnceDifferentiable object, and get reference to final U
    od, Uf = gen_opt_fun(Ut,Hd,Hc,t,n)
    # Optimization options
    #opt = Optim.Options(g_tol = 1E-9)
    # Run optimization
    res = optimize(od,vec(u_init.'),ConjugateGradient())
    #res = optimize(od,vec(u_init.'))
    # Reshape final control amplitude matrix
    uf = reshape(Optim.minimizer(res),m,n).'
    return uf, Uf, res
end

function opt_hadamard(n)
    Hd = σz
    Hc = [σx]
    t = 10.0
    #u_init = rand(n,1) .- 0.5
    u_init = zeros(n,1)
    Ut = Operator([1/√2 1/√2; 1/√2 -1/√2])
    grape(Ut,Hd,Hc,u_init,t,n)
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
    grape(Ut,Hd,Hc,u_init,t,n)
end

function plotgrape(tf,uf)
    figure()
    step(linspace(0,tf,size(uf,1)+1),[uf[1:1,:]; uf])
    legend(["Control $i" for i in 1:size(uf,2)])
    tight_layout()
    grid()
end
