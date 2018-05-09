# Schrodinger equation
function SchrodingerEvo(H₀::Operator)
    # Constant Hamiltonian term
    L₀ = -1im*data(H₀)
    return Liouvillian(dims(H₀),L₀)
end
function SchrodingerEvo(H₀::Operator, Hₙ::Tuple{Vararg{Tuple}})
    for Hᵢ in Hₙ; dimsmatch(H₀,first(Hᵢ)); end
    # Constant Hamiltonian term
    L₀ = -1im*data(H₀)
    # Time-dependent Hamiltonian terms
    Lₙ,fₙ,pₙ = unpack_operators(-1im,data,Hₙ)
    return Liouvillian(dims(H₀),L₀,Lₙ,fₙ,pₙ)
end
SchrodingerEvo(H₀::Operator,Hₙ::Vararg{Tuple}) = SchrodingerEvo(H₀,Hₙ)
SchrodingerEvo(H::Tuple{Operator,Vararg{Tuple}}) = SchrodingerEvo(first(H),tail(H))
SchrodingerEvo(::Any) = throw(ArgumentError("invalid Hamiltonian specification"))


# von Neumann equation
# TODO


# Lindblad form master equation
function LindbladEvo(H₀::Operator, Cₘ::Tuple{Vararg{Operator}})
    for Cᵢ in Cₘ; dimsmatch(H₀,Cᵢ); end
    I = data(qeye(dims(H₀)))
    # Constant Hamiltonian term
    L₀ = -1im*(I⊗data(H₀) - data(H₀).'⊗I)
    # Add constant collapse operator terms
    L₀ = add_constant_collapse(L₀,Cₘ,I)
    return Liouvillian(dims(H₀),L₀)
end
function LindbladEvo(H₀::Operator, Hₙ::Tuple{Vararg{Tuple}}, Cₘ::Tuple{Vararg{Operator}})
    for Hᵢ in Hₙ; dimsmatch(H₀,first(Hᵢ)); end
    for Cᵢ in Cₘ; dimsmatch(H₀,Cᵢ); end
    I = data(qeye(dims(H₀)))
    # Constant Hamiltonian term
    L₀ = -1im*(I⊗data(H₀) - transpose(data(H₀))⊗I)
    # Add constant collapse operator terms
    L₀ = add_constant_collapse(L₀,Cₘ,I)
    # Time-dependent Hamiltonian terms
    Lₙ,fₙ,pₙ = unpack_operators(-1im,data,Hₙ)
    return Liouvillian(dims(H₀),L₀,@.((I,)⊗Lₙ-transpose(Lₙ)⊗(I,)),fₙ,pₙ)
end
LindbladEvo(H₀::Operator,Cₘ::Vararg{Operator}) = LindbladEvo(H₀,Cₘ)
LindbladEvo(H::Tuple{Operator,Vararg{Tuple}},Cₘ::Vararg{Operator}) = LindbladEvo(first(H),tail(H),Cₘ)
LindbladEvo(H::Tuple{Operator,Vararg{Tuple}},Cₘ::Tuple{Vararg{Operator}}) = LindbladEvo(first(H),tail(H),Cₘ)
LindbladEvo(::Any) = throw(ArgumentError("invalid Hamiltonian specification"))


# Schrodinger evo propagator
SchrodingerProp(H₀::Operator, tspan) = SchrodingerProp(H₀,float(tspan[2]-tspan[1]))
function SchrodingerProp(H₀::Operator, Δt::Float64)
    # Constant Hamiltonian term
    U = expim(Hermitian(-full(H₀).*Δt))
    return Propagator(U,Δt,dims(H₀))
end
SchrodingerProp(H₀::Operator, Hₙ::Tuple, tspan, steps::Integer) = SchrodingerProp(H₀,(Hₙ,),tspan,steps)
function SchrodingerProp(H₀::Operator, Hₙ::Tuple{Vararg{Tuple}}, tspan, steps::Integer)
    for Hᵢ in Hₙ; dimsmatch(H₀,first(Hᵢ)); end
    # Sampling times and spacing dt
    t₁, t₂ = tspan; dt = (t₂-t₁)/steps
    # Unpack arguments
    Hn = unpack_operators(1,full,Hₙ)
    # Constant Hamiltonian term
    U₀ = expim(Hermitian(-dt*full(H₀)))
    # Multiply-in sampled time-dependent terms interleaved with constant term
    U = eye(U₀); H = Matrix{Base.promote_eltype(Hn[1]...)}(size(H₀)...)
    A = similar(U)
    for i = 1:steps
        A_mul_B!(A,U₀,U); copy!(U,A)
        step_hamiltonian!(H,Hn,(t₁,dt,i))
        expim!(A,Hermitian(H))
        U = A*U
    end
    return Propagator(U,float(t₂-t₁),dims(H₀))
end
SchrodingerProp(H) = throw(ArgumentError("invalid Propagator specification"))


# Lindblad evo propagator
LindbladProp(H₀::Operator, Cₘ::Operator, tspan) = LindbladProp(H₀,(Cₘ,),tspan)
LindbladProp(H₀::Operator, Cₘ::Tuple{Vararg{Operator}}, tspan) = LindbladProp(H₀,Cₘ,tspan[2]-tspan[1])
function LindbladProp(H₀::Operator, Cₘ::Tuple{Vararg{Operator}}, Δt::Float64)
    for Cᵢ in Cₘ; dimsmatch(H₀,Cᵢ); end
    I = eye(prod(dims(H₀)))
    # Constant Hamiltonian term
    L₀ = -1im.*(I⊗full(H₀) .- transpose(full(H₀))⊗I)
    # Constant collapse operator terms
    L₀ = add_constant_collapse(L₀,Cₘ,I)
    # Build constant propagator
    U = LinAlg.expm!(L₀*Δt)
    return Propagator(U,Δt,dims(H₀))
end
LindbladProp(H₀::Operator, Hₙ::Tuple, Cₘ::Operator, tspan, steps::Integer) = LindbladProp(H₀,(Hₙ,),(Cₘ,),tspan,steps)
LindbladProp(H₀::Operator, Hₙ::Tuple, Cₘ::Tuple{Vararg{Operator}}, tspan, steps::Integer) = LindbladProp(H₀,(Hₙ,),Cₘ,tspan,steps)
LindbladProp(H₀::Operator, Hₙ::Tuple{Vararg{Tuple}}, Cₘ::Operator, tspan, steps::Integer) = LindbladProp(H₀,Hₙ,(Cₘ,),tspan,steps)
function LindbladProp(H₀::Operator, Hₙ::Tuple{Vararg{Tuple}}, Cₘ::Tuple{Vararg{Operator}}, tspan, steps::Integer)
    for Hᵢ in Hₙ; dimsmatch(H₀,first(Hᵢ)); end
    for Cᵢ in Cₘ; dimsmatch(H₀,Cᵢ); end
    I = eye(prod(dims(H₀)))
    # Sampling times and spacing dt
    t₁, t₂ = tspan; dt = (t₂-t₁)/steps
    # Unpack arguments
    Hn = unpack_operators(1,full,Hₙ)
    # Constant Hamiltonian term
    L₀ = -1im.*(I⊗full(H₀) .- transpose(full(H₀))⊗I)
    # Constant collapse operator terms
    L₀ = add_constant_collapse(L₀,Cₘ,I)
    # Build constant propagator part
    U₀ = LinAlg.expm!(L₀.*dt)
    # Multiply-in sampled time-dependent terms interleaved with constant term
    U = eye(U₀); H = Matrix{Base.promote_eltype(Hn[1]...)}(size(H₀)...)
    A = similar(U); B = Matrix{Complex128}(size(H₀)...)
    for i = 1:steps
        A_mul_B!(A,U₀,U); copy!(U,A)
        step_hamiltonian!(H,Hn,(t₁,dt,i))
        expim!(B,Hermitian(H))
        invB = LinAlg.inv!(lufact(B))
        At_mul_B!(A,invB⊗I,U)
        I_kron_A_mul_B!(U,B,A) # A_mul_B!(U,Id⊗A,C)
    end
    return Propagator(U,float(t₂-t₁),dims(H₀))
end
LindbladProp(H) = throw(ArgumentError("invalid Propagator specification"))


# Utils
function unpack_operators(x::Number,f::Function,t::Tuple{Tuple{Operator,Function,Array},Vararg{Tuple}})
    a, b = first(t), unpack_operators(tail(t))
    return (x*f(a[1]), b[1]...), (a[2], b[2]...), (a[3], b[3]...)
end
function unpack_operators(x::Number,f::Function,t::Tuple{Tuple{Operator,Function},Vararg{Tuple}})
    a, b = first(t), unpack_operators(tail(t))
    return (x*f(a[1]), b[1]...), (a[2], b[2]...), ([], b[3]...)
end
unpack_operators(t::Tuple{}) = (), (), ()

function add_constant_collapse(L₀,Cₘ,I)
    for Cᵢ in Cₘ
        C = data(Cᵢ); CdC = C'*C
        L₀ += conj(C)⊗C - 0.5*(I⊗CdC + CdC.'⊗I)
    end
    return L₀
end

function step_hamiltonian!(H,Hn,tspec)
    Hₙ,fₙ,pₙ = Hn
    t₁, dt, i = tspec; tᵢ = t₁ + (i-1)*dt
    fill!(H,0) # zero out matrix
    for j = 1:length(Hₙ) # sample function at 3 points
        Hᵢ,fᵢ,pᵢ = Hₙ[j],fₙ[j],pₙ[j]
        H .+= -dt*(fᵢ(tᵢ,pᵢ) + fᵢ(tᵢ+0.5dt,pᵢ) + fᵢ(tᵢ+dt,pᵢ))/3 .* Hᵢ
    end
end
