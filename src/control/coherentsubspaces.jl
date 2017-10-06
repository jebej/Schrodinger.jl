immutable CoherentSubspaces{T,D} <: ObjectiveFunction
    dims::SDims{D}
    δt::Float64 # timestep (fixed)
    Ut::Matrix{Complex128} # target unitary
    s::IntSet # coherent subspaces
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

CoherentSubspaces(dims,δt,Ut,s,Hd,Hc,u_last,U,X,P,D,V,H,A,Jkj,cisDj) = CoherentSubspaces{eltype(Hd),length(dims)}(dims,δt,Ut,s,Hd,Hc,u_last,U,X,P,D,V,H,A,Jkj,cisDj)

function CoherentSubspaces(Ut::Operator,s::IntCol,Hd::Operator,Hc::Vector{<:Operator},t::Real,n::Integer)
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
    return CoherentSubspaces(dims(Hd),t/n,Ut_d,s,Hd_d,Hc_d,u_last,U,X,P,D,V,H,A,Jkj,cisDj)
end

function objective(O::CoherentSubspaces,u)
    N = size(O.Ut,1); Uf = O.X[end]
    # Calculate forward propagators
    calc_fprops!(O,u)
    # Calculate coherent subspaces norm infidelity:
    # fₑ = 1 - |⟨Ut,Uf⟩|cs/N
    return 1 - inner_cs(O.Ut,Uf,O.s)/N
end

function gradient!(O::CoherentSubspaces,fp,u)
    n = length(O.U); m = length(O.Hc); N = size(O.Ut,1); Uf = O.X[end]
    # Calculate forward and backward propagators
    calc_fprops!(O,u)
    calc_bprops!(O)
    # Calculate exact derivative of coherent subspaces norm infidelity
    x,y = _inner_cs_1(Uf,O.Ut,O.s)
    for j = 1:n
        O.cisDj .= cis.(O.D[j])
        for k = 1:m
            Jmat!(O.Jkj,O.Hc[k],O.cisDj,O.D[j],O.V[j],O.δt,O.A)
            j==1 ? copy!(O.A,O.Jkj) : A_mul_B!(O.A,O.Jkj,O.X[j-1])
            fp[(j-1)*m+k] = -inner_cs_grad(O.P[j],O.A,x,y,O.s)/N
        end
    end
    return fp
end
