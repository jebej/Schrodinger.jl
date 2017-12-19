immutable NormPSU{T,D} <: ObjectiveFunction
    dims::SDims{D}
    δt::Float64 # timestep (fixed)
    Ut::Matrix{Complex128} # target unitary
    Hd::Matrix{T} # drift Hamiltonian
    Hc::Vector{Matrix{T}} # control Hamiltonian(s)
    u_last::Vector{Float64} # last control amplitudes, ordered by control, and then time
    # Propagator and eigendecompotion storage
    U::Vector{Matrix{Complex128}} # individual propagators
    X::Vector{Matrix{Complex128}} # forward cumulative propagators
    P::Vector{Matrix{Complex128}} # backwards cumulative propagators
    D::Vector{Vector{Float64}} # eigenvalues associated to each propagator
    V::Vector{Matrix{T}} # eigenvector matrix associated to each propagator
    # Workspace variables
    H::Matrix{T} # temporary storage for full Hamiltonian
    A::Matrix{Complex128} # temporary
    Jkj::Matrix{Complex128} # temporary
    cisDj::Vector{Complex128} # temporary
end

NormPSU(dims,δt,Ut,Hd,Hc,u_last,U,X,P,D,V,H,A,Jkj,cisDj) = NormPSU{eltype(H),length(dims)}(dims,δt,Ut,Hd,Hc,u_last,U,X,P,D,V,H,A,Jkj,cisDj)

function NormPSU(Ut::Operator,Hd::Operator,Hc::Vector{<:Operator},t::Real,n::Integer)
    # TODO: Prob should do some error checking
    N = prod(dims(Hd))
    m = length(Hc)
    # Make sure we pass dense operators
    Ut_d = full(Ut)
    Hd_d = full(Hd)
    Hc_d = full.(Hc)
    # Storage for last control ampitudes, NaN for first run
    u_last = fill(NaN64,m*n)
    # Generate cache for various objects
    H = promote_type(typeof(Hd_d),eltype(Hc_d))(N,N)
    A = Matrix{Complex128}(N,N)
    U = [Matrix{Complex128}(N,N) for i=1:n]
    X = [Matrix{Complex128}(N,N) for i=1:n]
    P = [Matrix{Complex128}(N,N) for i=1:n]
    # For the exact derivative we need to store eigenvectors and eigenvalues
    D = [Vector{Float64}(N) for i=1:n] # eigenvalues are always real
    V = [similar(H) for i=1:n]
    # More temporary storage
    Jkj = Matrix{Complex128}(N,N)
    cisDj = Vector{Complex128}(N)
    return NormPSU(dims(Hd),t/n,Ut_d,Hd_d,Hc_d,u_last,U,X,P,D,V,H,A,Jkj,cisDj)
end

function objective(O::NormPSU,u)
    N = size(O.Ut,1); Uf = O.X[end]
    # Calculate forward propagators
    calc_fprops!(O,u)
    # Calculate PSU (projective special unitary) norm infidelity:
    # fₑ = 1 - |⟨Ut,Uf⟩|/N
    return 1 - abs(inner(O.Ut,Uf))/N
end

function gradient!(O::NormPSU,fp,u)
    n = length(O.U); m = length(O.Hc); N = size(O.Ut,1); Uf = O.X[end]
    # Calculate forward and backward propagators
    calc_fprops!(O,u)
    calc_bprops!(O)
    # Calculate exact derivative of fidelity error function:
    # ∂fₑ/∂uₖⱼ = -1/2*⟨Pⱼ,∂Uⱼ/∂uₖⱼ*Xⱼ₋₁⟩*cis(∠⟨Xⱼ,Pⱼ⟩)/N + c.c.
    #          = -Re(⟨Pⱼ,Jₖⱼ*Xⱼ₋₁⟩*cis(∠⟨Uf,Ut⟩))/N where Jₖⱼ = ∂Uⱼ/∂uₖⱼ
    a = normalize(inner(Uf,O.Ut))
    for j = 1:n
        O.cisDj .= cis.(O.D[j])
        for k = 1:m
            Jmat!(O.Jkj,O.Hc[k],O.cisDj,O.D[j],O.V[j],O.δt,O.A)
            j==1 ? copy!(O.A,O.Jkj) : A_mul_B!(O.A,O.Jkj,O.X[j-1])
            fp[(k-1)*n+j] = -real(inner(O.P[j],O.A)*a)/N
        end
    end
    return fp
end
