struct CoherentSubspaces{M,S<:Union{Matrix{ComplexF64},NTuple{M,Matrix{ComplexF64}}},T,D} <: ObjectiveFunction
    dims::Dims{D}
    δt::Float64 # timestep (fixed)
    Ut::S # target unitary
    s::BitSet # coherent subspaces
    Ne::Int # effective target Operator dimension (rank of target operator)
    Hd::Matrix{T} # drift Hamiltonian
    Hc::Vector{Matrix{T}} # control Hamiltonian(s)
    u_last::Vector{Float64} # last control amplitudes, ordered by control, and then time
    # Propagator and eigendecompotion storage
    U::Vector{Matrix{ComplexF64}} # individual propagators
    X::Vector{Matrix{ComplexF64}} # forward cumulative propagators
    P::Vector{S} # backwards cumulative propagators
    D::Vector{Vector{Float64}} # eigenvalues associated to each propagator
    V::Vector{Matrix{T}} # eigenvector matrix associated to each propagator
    # Workspace variables
    H::Hermitian{T,Matrix{T}} # temporary storage for full Hamiltonian
    A::Matrix{ComplexF64} # temporary
    Jkj::Matrix{ComplexF64} # temporary
    cisDj::Vector{ComplexF64} # temporary
end

CoherentSubspaces(dims,δt,Ut::NTuple{M,<:Matrix},s,Ne,Hd,Hc,u_last,U,X,P,D,V,H,A,Jkj,cisDj) where M = CoherentSubspaces{M,typeof(Ut),eltype(H),length(dims)}(dims,δt,Ut,s,Ne,Hd,Hc,u_last,U,X,P,D,V,H,A,Jkj,cisDj)
CoherentSubspaces(dims,δt,Ut::Matrix,s,Ne,Hd,Hc,u_last,U,X,P,D,V,H,A,Jkj,cisDj) = CoherentSubspaces{1,typeof(Ut),eltype(H),length(dims)}(dims,δt,Ut,s,Ne,Hd,Hc,u_last,U,X,P,D,V,H,A,Jkj,cisDj)

function CoherentSubspaces(Ut::NTuple{M,<:Operator},s::IntCol,Hd::Operator,Hc::Vector{<:Operator},t::Real,n::Integer,Ne=prod(dims(Hd))) where M
    # TODO: Prob should do some error checking on inputs
    N = prod(dims(Hd))
    m = length(Hc)
    # Make sure we pass dense operators
    Ut_d = complex.(full.(Ut))
    Hd_d = full(Hd)
    Hc_d = full.(Hc)
    # Storage for last control ampitudes, NaN for first run
    u_last = fill(NaN64,m*n)
    # Generate cache for various objects
    H = Hermitian(zeros(promote_eltype(typeof(Hd_d),eltype(Hc_d)),N,N))
    A = Matrix{ComplexF64}(N,N)
    U = [Matrix{ComplexF64}(N,N) for i=1:n]
    X = [Matrix{ComplexF64}(N,N) for i=1:n]
    P = [similar.(Ut_d) for i=1:n]
    for i=1:M; copy!(P[end][i],Ut_d[i]); end
    # For the exact derivative we need to store eigenvectors and eigenvalues
    D = [Vector{Float64}(N) for i=1:n] # eigenvalues are always real
    V = [similar(H.data) for i=1:n]
    # More temporary storage
    Jkj = Matrix{ComplexF64}(N,N)
    cisDj = Vector{ComplexF64}(N)
    return CoherentSubspaces(dims(Hd),t/n,Ut_d,s,Ne,Hd_d,Hc_d,u_last,U,X,P,D,V,H,A,Jkj,cisDj)
end

function CoherentSubspaces(Ut::Operator,s::IntCol,Hd::Operator,Hc::Vector{<:Operator},t::Real,n::Integer,Ne=prod(dims(Hd)))
    # TODO: Prob should do some error checking on inputs
    N = prod(dims(Hd))
    m = length(Hc)
    # Make sure we pass dense operators
    Ut_d = complex.(full.(Ut))
    Hd_d = full(Hd)
    Hc_d = full.(Hc)
    # Storage for last control ampitudes, NaN for first run
    u_last = fill(NaN64,m*n)
    # Generate cache for various objects
    H = Hermitian(zeros(promote_eltype(typeof(Hd_d),eltype(Hc_d)),N,N))
    A = Matrix{ComplexF64}(N,N)
    U = [Matrix{ComplexF64}(N,N) for i=1:n]
    X = [Matrix{ComplexF64}(N,N) for i=1:n]
    P = [Matrix{ComplexF64}(N,N) for i=1:n]
    copy!(P[end],Ut_d)
    # For the exact derivative we need to store eigenvectors and eigenvalues
    D = [Vector{Float64}(N) for i=1:n] # eigenvalues are always real
    V = [similar(H.data) for i=1:n]
    # More temporary storage
    Jkj = Matrix{ComplexF64}(N,N)
    cisDj = Vector{ComplexF64}(N)
    return CoherentSubspaces(dims(Hd),t/n,Ut_d,s,Ne,Hd_d,Hc_d,u_last,U,X,P,D,V,H,A,Jkj,cisDj)
end

function objective(O::CoherentSubspaces,u)
    Uf = O.X[end]
    # Calculate forward propagators
    calc_fprops!(O,u)
    # Calculate coherent subspaces norm infidelity:
    # fₑ = 1 - |⟨Ut,Uf⟩|cs/N
    return 1 - inner_cs(O.Ut,Uf,O.s)/O.Ne
end

function gradient!(O::CoherentSubspaces,fp,u)
    Uf = O.X[end]
    n = length(O.U); m = length(O.Hc)
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
            fp[(k-1)*n+j] = -inner_cs_grad(O.P[j],O.A,x,y,O.s)/O.Ne
        end
    end
    return fp
end
