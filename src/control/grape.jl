using Schrodinger, Optim

function calc_fprops!(U,X,Hd,Hc,u,δt,H,u_last)
    # If the control amplitudes u did not change, return
    u == u_last && return nothing
    # Otherwise calculate each individual propagator
    n = length(U); m = length(Hc)
    for j = 1:n
        copy!(H,Hd)
        for k = 1:m
            H .+= u[(j-1)*m+k].*Hc[k]
        end
        Schrodinger.expim!(U[j],Hermitian(scale!(H,-δt)))
    end
    # Calculate cumulative product
    copy!(X[1],U[1])
    for i=2:n
        A_mul_B!(X[i],U[i],X[i-1])
    end
    # Save control amplitudes
    copy!(u_last,u)
    return nothing
end

function calc_bprops!(P,U,Ut)
    copy!(P[end],Ut)
    for i=length(P)-1:-1:1
        Ac_mul_B!(P[i],U[i+1],P[i+1])
    end
    return nothing
end

function hsinner!(U1::AbstractArray,U2::AbstractArray,C::AbstractArray)
    Ac_mul_B!(C,U1,U2)
    return trace(C)
end

function fidelity!(Ut,Hd,Hc,u,δt,A,U,X,u_last)
    calc_fprops!(U,X,Hd,Hc,u,δt,A,u_last)
    # Calculate fidelity
    return 1-abs2(hsinner!(Ut,X[end],A))/size(Ut,1)^2
end

function fidelityprime!(fp,Ut,Hd,Hc,u,δt,A,U,X,P,u_last)
    calc_fprops!(U,X,Hd,Hc,u,δt,A,u_last)
    calc_bprops!(P,U,Ut)
    n = length(U); m = length(Hc)
    # Calculate derivative approximation
    for j = 1:n, k=1:m
        fp[(j-1)*m+k] = real(hsinner!(P[j],scale!(Hc[k]*X[j],1im*δt),A)*hsinner!(X[j],P[j],A))
    end
    scale!(fp,2/size(Ut,1)^2)
    return fp
end

function gen_opt_fun(Ut::Operator,Hd::Operator,Hc::Vector{<:Operator},t::Real,n::Integer)
    N = prod(dims(Hd))
    m = length(Hc)
    # Generate cache for various objects
    A = Matrix{Complex128}(N,N)
    U = [Matrix{Complex128}(N,N) for i=1:n]
    X = deepcopy(U)
    P = deepcopy(U)
    u_last = Vector{Float64}(m*n)
    # Make sure we pass dense operators
    Ut_d = full(Ut)
    Hd_d = full(Hd)
    Hc_d = [full(H) for H in Hc]
    # Create optimization function
    f = (u) -> fidelity!(Ut_d,Hd_d,Hc_d,u,t/n,A,U,X,u_last)
    g! = (fp,u) -> fidelityprime!(fp,Ut_d,Hd_d,Hc_d,u,t/n,A,U,X,P,u_last)
    # Return a OnceDifferentiable object with appropriate seed
    return OnceDifferentiable(f,g!,zeros(m*n)),X[end]
end


function grape(Ut::Operator,Hd::Operator,Hc::Vector{<:Operator},u_init,t::Real,n::Integer)
    m = length(Hc)
    (n,m) == size(u_init) || throw(ArgumentError("control amplitude matrix not consistent with number of timesteps or number of control hamiltonians"))
    # Generate OnceDifferentiable object, and get reference to final U
    od,Uf = gen_opt_fun(Ut,Hd,Hc,t,n)
    # Run optimization
    res = optimize(od,vec(u_init.'),ConjugateGradient())
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
    Hd = 0.5 * complex(σx⊗σx + σy⊗σy +σz⊗σz)
    Hc = [0.5*complex(σx⊗σ0), 0.5*σy⊗σ0, 0.5*complex(σ0⊗σx), 0.5*σ0⊗σy]
    t = 6.0
    #u_init = zeros(n,4)
    u_init = rand(n,4) .- 0.5
    Ut = Operator([0.5  0.5    0.5  0.5
                   0.5  0.5im -0.5 -0.5im
                   0.5 -0.5    0.5 -0.5
                   0.5 -0.5im -0.5  0.5im],
                  (2,2))
    grape(Ut,Hd,Hc,u_init,t,n)
end
