
function infidelity!(u,δt,Ut,Hd,Hc,H,U,X,D,V,u_last)
    # Calculate forward propagators
    calc_fprops!(U,X,D,V,u,δt,Hd,Hc,H,u_last)
    Uf = X[end]; N² = size(Ut,1)^2
    # Calculate infidelity:
    # fₑ = 1 - Φ where Φ = |⟨Ut,Uf⟩|²/N²
    return 1 - abs2(inner(Ut,Uf))/N²
end

function infidelityprime!(fp,u,δt,Ut,Hd,Hc,H,A,U,X,D,V,P,u_last)
    # Calculate forward and backward propagators
    calc_fprops!(U,X,D,V,u,δt,Hd,Hc,H,u_last)
    calc_bprops!(P,U,Ut)
    n = length(U); m = length(Hc); N² = size(Ut,1)^2
    # Calculate exact derivative of the infidelity function:
    # ∂Φ/∂uₖⱼ  = ⟨Pⱼ,∂Uⱼ/∂uₖⱼ*Xⱼ₋₁⟩⟨Xⱼ,Pⱼ⟩/N² + c.c.
    # ∂fₑ/∂uₖⱼ = -∂Φ/∂uₖⱼ
    #          = -2*Re(⟨Pⱼ,Jₖⱼ*Xⱼ₋₁⟩⟨Uf,Ut⟩)/N² where Jₖⱼ = ∂Uⱼ/∂uₖⱼ
    Jkj = similar(A)
    cisDj = Vector{ComplexF64}(length(D[1]))
    a = inner(X[end],Ut)
    for j = 1:n
        cisDj .= cis.(D[j])
        for k = 1:m
            Jmat!(Jkj,Hc[k],cisDj,D[j],V[j],δt,A)
            j==1 ? copy!(A,Jkj) : mul!(A,Jkj,X[j-1])
            fp[(k-1)*n+j] = -2*real(inner(P[j],A)*a)/N²
        end
    end
    return fp
end

function infidelityprimeapprox!(fp,u,δt,Ut,Hd,Hc,H,A,U,X,D,V,P,u_last)
    # Calculate forward and backward propagators
    calc_fprops!(U,X,D,V,u,δt,Hd,Hc,H,u_last)
    calc_bprops!(P,U,Ut)
    n = length(U); m = length(Hc); N² = size(Ut,1)^2
    # Calculate approximate derivative of fidelity error function:
    # ∂Φ/∂uₖⱼ  = ⟨Pⱼ,∂Uⱼ/∂uₖⱼ*Xⱼ₋₁⟩⟨Xⱼ,Pⱼ⟩/N² + c.c.
    # ∂fₑ/∂uₖⱼ = -∂Φ/∂uₖⱼ
    #          ≈ 2*Re(⟨Pⱼ,iδt*Hₖ*Xⱼ⟩⟨Uf,Ut⟩)/N²
    a = inner(X[end],Ut)
    for j = 1:n, k = 1:m
        mul!(A,Hc[k],X[j])
        fp[(k-1)*n+j] = -2δt*imag(inner(P[j],A)*a)/N²
    end
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
    A = Matrix{ComplexF64}(N,N)
    U = [similar(A) for i=1:n]
    X = deepcopy(U)
    P = deepcopy(U)
    # For the exact derivative we need to store eigenvectors and eigenvalues
    V = [similar(H) for i=1:n]
    D = [Vector{Float64}(N) for i=1:n] # eigenvalues are always real
    # Storage for last control ampitudes, NaN for first run
    u_last = fill(NaN64,m*n)
    # Create optimization function object
    f = (u) -> infidelity!(u,t/n,Ut_d,Hd_d,Hc_d,H,U,X,D,V,u_last)
    g! = (fp,u) -> infidelityprime!(fp,u,t/n,Ut_d,Hd_d,Hc_d,H,A,U,X,D,V,P,u_last)
    #g! = (fp,u) -> infidelityprimeapprox!(fp,u,t/n,Ut_d,Hd_d,Hc_d,H,A,U,X,D,V,P,u_last)
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
    res = optimize(od,vec(transpose(u_init)),ConjugateGradient())
    #res = optimize(od,vec(transpose(u_init)))
    # Reshape final control amplitude matrix
    uf = transpose(reshape(Optim.minimizer(res),m,n))
    return uf, Uf, res
end

function opt_2Q_QFT1(n)
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
