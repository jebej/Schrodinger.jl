# Schrodinger equation
function SchrodingerEvo(H₀::Operator)
    D = length(dims(H₀))
    # Constant Hamiltonian term
    L₀ = -1.0im.*data(H₀)
    return Liouvillian{0,D}(dims(H₀),L₀)
end
function SchrodingerEvo(H₀::Operator, Hₙ::Vector)
    N = length(Hₙ); D = length(dims(H₀));
    for i=1:N; dimsmatch(H₀,Hₙ[i][1]); end
    # Constant Hamiltonian term
    L₀ = -1.0im.*data(H₀)
    # Time-dependent Hamiltonian terms
    Lₙ = ((-1.0im.*data(H[1]) for H in Hₙ)...)
    fₙ = ((H[2] for H in Hₙ)...)
    pₙ = ((length(H)==3?H[3]:[] for H in Hₙ)...)
    return Liouvillian{N,D}(dims(H₀),L₀,Lₙ,fₙ,pₙ)
end
function SchrodingerEvo(H::Vector)
    H₀ = H[1]::Operator
    Hₙ = H[2:end]::Vector
    return SchrodingerEvo(H₀,Hₙ)
end
SchrodingerEvo(H) = throw(ArgumentError("invalid Hamiltonian specification"))


# von Neumann equation
# TODO


# Lindblad form master equation
function LindbladEvo(H₀::Operator, Cₘ::Vector)
    M = length(Cₘ); D = length(dims(H₀));
    for i=1:M; dimsmatch(H₀,Cₘ[i]); end
    Id = data(qeye(prod(dims(H₀))))
    # Constant Hamiltonian term
    L₀ = -1.0im.*(Id⊗data(H₀) - data(H₀).'⊗Id)
    # Constant collapse operator terms
    for i = 1:M
        C = data(Cₘ[i])
        CdC = C'*C
        L₀ .+= conj(C)⊗C - 0.5.*(Id⊗CdC + CdC.'⊗Id)
    end
    return Liouvillian{0,D}(dims(H₀),L₀)
end
function LindbladEvo(H₀::Operator, Hₙ::Vector, Cₘ::Vector)
    N = length(Hₙ); M = length(Cₘ); D = length(dims(H₀));
    for i=1:N; dimsmatch(H₀,Hₙ[i][1]); end
    for i=1:M; dimsmatch(H₀,Cₘ[i]); end
    Id = data(qeye(prod(dims(H₀))))
    # Constant Hamiltonian term
    L₀ = -1.0im.*(Id⊗data(H₀) - data(H₀).'⊗Id)
    # Constant collapse operator terms
    for i = 1:M
        C = data(Cₘ[i])
        CdC = C'*C
        L₀ .+= conj(C)⊗C - 0.5.*(Id⊗CdC + CdC.'⊗Id)
    end
    # Time-dependent Hamiltonian terms
    Lₙ = ((-1.0im.*(Id⊗data(H[1]) - data(H[1]).'⊗Id) for H in Hₙ)...)
    fₙ = ((H[2] for H in Hₙ)...)
    pₙ = ((length(H)==3?H[3]:[] for H in Hₙ)...)
    return Liouvillian{N,D}(dims(H₀),L₀,Lₙ,fₙ,pₙ)
end
function LindbladEvo(H::Vector, Cₘ::Vector)
    H₀ = H[1]::Operator
    Hₙ = H[2:end]::Vector
    return LindbladEvo(H₀,Hₙ,Cₘ)
end
LindbladEvo(H,vararg...) = throw(ArgumentError("invalid Hamiltonian specification"))


# Schrodinger evo propagator
function SchrodingerProp(H₀::Operator, Δt::Float64)
    D = length(dims(H₀))
    # Constant Hamiltonian term
    U = expim(Hermitian(-full(H₀).*Δt))
    return Propagator{D}(dims(H₀),U,Δt)
end
function SchrodingerProp(H₀::Operator, Hₙ::Vector, Δt::Float64, n::Int)
    N = length(Hₙ); D = length(dims(H₀));
    for i=1:N; dimsmatch(H₀,Hₙ[i][1]); end
    # Sampling times and spacing dt
    ts = linspace(0,Δt,n+1)[2:end]; dt = ts[2]-ts[1]
    # Constant Hamiltonian term
    U₀ = expim(Hermitian(-full(H₀).*dt))
    # Multiply-in sampled time-dependent terms interleaved with constant term
    U = eye(Complex128,prod(dims(H₀)))
    for t in ts
        U = U₀*U
        H = zeros(U)
        for i in 1:N
            f = Hₙ[i][2]
            p = length(Hₙ[i])==3?Hₙ[i][3]:[]
            H .+= full(Hₙ[i][1]).*(f(t,p)*dt)
        end
        U = expim(Hermitian(-H))*U
    end
    return Propagator{D}(dims(H₀),U,Δt)
end
function SchrodingerProp(H::Vector, Δt::Float64, n::Int)
    H₀ = H[1]::Operator
    Hₙ = H[2:end]::Vector
    return SchrodingerProp(H₀,Hₙ,Δt,n)
end
SchrodingerProp(H) = throw(ArgumentError("invalid Propagator specification"))


# Lindblad evo propagator
function LindbladProp(H₀::Operator, Cₘ::Vector, Δt::Float64)
    M = length(Cₘ); D = length(dims(H₀));
    for i=1:M; dimsmatch(H₀,Cₘ[i]); end
    Id = data(qeye(prod(dims(H₀))))
    # Constant Hamiltonian term
    L₀ = -1.0im.*(Id⊗data(H₀) - data(H₀).'⊗Id)
    # Constant collapse operator terms
    for i = 1:M
        C = data(Cₘ[i])
        CdC = C'*C
        L₀ += conj(C)⊗C - 0.5.*(Id⊗CdC + CdC.'⊗Id)
    end
    # Build constant propagator
    U = expm(full(L₀).*Δt)
    return Propagator{D}(dims(H₀),U,Δt)
end
function LindbladProp(H₀::Operator, Hₙ::Vector, Cₘ::Vector, Δt::Float64, n::Int)
    N = length(Hₙ); M = length(Cₘ); D = length(dims(H₀));
    for i=1:N; dimsmatch(H₀,Hₙ[i][1]); end
    for i=1:M; dimsmatch(H₀,Cₘ[i]); end
    # Sampling times and spacing dt
    ts = linspace(0,Δt,n+1)[2:end]; dt = ts[2]-ts[1]
    Id = data(qeye(prod(dims(H₀))))
    # Constant Hamiltonian term
    L₀ = -1.0im.*(Id⊗data(H₀) - data(H₀).'⊗Id)
    # Constant collapse operator terms
    for i = 1:M
        C = data(Cₘ[i])
        CdC = C'*C
        L₀ += conj(C)⊗C - 0.5.*(Id⊗CdC + CdC.'⊗Id)
    end
    # Build constant propagator part
    U₀ = LinAlg.expm!(full(L₀).*dt)
    # Multiply-in sampled time-dependent terms interleaved with constant term
    Lₙ = ((full(H[1]) for H in Hₙ)...)
    fₙ = ((H[2] for H in Hₙ)...)
    pₙ = ((length(H)==3?H[3]:[] for H in Hₙ)...)
    U = eye(Complex128,prod(dims(H₀))^2)
    L = Matrix{Complex128}(prod(dims(H₀)),prod(dims(H₀)))
    C = similar(U)
    for t in ts
        A_mul_B!(C,U₀,U); copy!(U,C) # U = U₀*U, slow for large matrices
        fill!(L,0)
        for i in 1:N
            L .-= Lₙ[i].*(fₙ[i](t,pₙ[i])*dt)
        end
        A = expim!(Hermitian(L))
        invA = LinAlg.inv!(lufact(A))
        At_mul_B!(C,invA⊗Id,U)
        I_kron_A_mul_B!(U,A,C) # A_mul_B!(U,Id⊗A,C)
    end
    return Propagator{D}(dims(H₀),U,Δt)
end
function LindbladProp(H::Vector, Cₘ::Vector, Δt::Float64, n::Int)
    H₀ = H[1]::Operator
    Hₙ = H[2:end]::Vector
    return LindbladProp(H₀,Hₙ,Cₘ,Δt,n)
end
LindbladProp(H) = throw(ArgumentError("invalid Propagator specification"))
