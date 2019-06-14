# Schrodinger equation
function SchrodingerEvo(H₀::Operator)
    # Constant Hamiltonian term
    L₀ = -1im*data(H₀)
    return Liouvillian(dims(H₀),L₀)
end
function SchrodingerEvo(H₀::Operator, Hₙ::Tuple{Vararg{Tuple}})
    dimsmatch(H₀,first.(Hₙ))
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
    dimsmatch(H₀,Cₘ)
    Id = data(qeye(dims(H₀)))
    # Constant Hamiltonian term
    L₀ = -1im*(Id⊗data(H₀) - transpose(data(H₀))⊗Id)
    # Add constant collapse operator terms
    L₀ += sum_collapse(Cₘ,Id,1)
    return Liouvillian(dims(H₀),L₀)
end
function LindbladEvo(H₀::Operator, Hₙ::Tuple{Vararg{Tuple}}, Cₘ::Tuple{Vararg{Operator}})
    dimsmatch(H₀,first.(Hₙ)); dimsmatch(H₀,Cₘ)
    Id = data(qeye(dims(H₀)))
    # Constant Hamiltonian term
    L₀ = -1im*(Id⊗data(H₀) - transpose(data(H₀))⊗Id)
    # Add constant collapse operator terms
    L₀ += sum_collapse(Cₘ,Id,1)
    # Time-dependent Hamiltonian terms
    Lₙ,fₙ,pₙ = unpack_operators(-1im,data,Hₙ)
    return Liouvillian(dims(H₀),L₀,@.((Id,)⊗Lₙ-transpose(Lₙ)⊗(Id,)),fₙ,pₙ)
end
LindbladEvo(H₀::Operator,Cₘ::Vararg{Operator}) = LindbladEvo(H₀,Cₘ)
LindbladEvo(H::Tuple{Operator,Vararg{Tuple}},Cₘ::Vararg{Operator}) = LindbladEvo(first(H),tail(H),Cₘ)
LindbladEvo(H::Tuple{Operator,Vararg{Tuple}},Cₘ::Tuple{Vararg{Operator}}) = LindbladEvo(first(H),tail(H),Cₘ)
LindbladEvo(::Any) = throw(ArgumentError("invalid Hamiltonian specification"))


# Schrodinger evo propagator
SchrodingerProp(H₀::Operator, tspan) = SchrodingerProp(H₀,float(tspan[2]-tspan[1]))
function SchrodingerProp(H₀::Operator, Δt::Real)
    # Constant Hamiltonian term
    U = expim(Hermitian(-full(H₀).*Δt))
    return Propagator(U,Δt,dims(H₀))
end
SchrodingerProp(H₀::Operator, Hₙ::Tuple, tspan, steps::Integer) = SchrodingerProp(H₀,(Hₙ,),tspan,steps)
function SchrodingerProp(H₀::Operator, Hₙ::Tuple{Vararg{Tuple}}, tspan, steps::Integer)
    dimsmatch(H₀,first.(Hₙ))
    F = real(eltype(H₀))
    # Sampling times and spacing dt
    t₁, t₂ = tspan; dt = (t₂-t₁)/steps
    # Unpack constant and time dep operators
    H0, Hn = full(H₀), unpack_operators(1,full,Hₙ)
    # Multiply sampled propagators together to generate total evolution
    U = Matrix{Complex{F}}(I, size(H0))
    H = Hermitian(zeros(compute_H_type(H0,Hn),size(H0)...))
    A = similar(U); B = similar(U); C = similar(H.data); D = similar(U)
    Λ = Vector{F}(undef,size(H,1))
    for i = 1:steps
        step_hamiltonian!(H.data,H0,Hn,(t₁,dt,i)) # calc H for this time step
        expim!(A,H,Λ,C,D) # A = exp(-1im*H*dt)
        mul!(B,A,U) # U = A*U
        U,B = B,U # swappitty swap for the next step
    end
    return Propagator(U,float(t₂-t₁),dims(H₀))
end
SchrodingerProp(H) = throw(ArgumentError("invalid Propagator specification"))


# Lindblad evo propagator
LindbladProp(H₀::Operator, Cₘ::Operator, tspan) = LindbladProp(H₀,(Cₘ,),tspan)
LindbladProp(H₀::Operator, Cₘ::Tuple{Vararg{Operator}}, tspan) = LindbladProp(H₀,Cₘ,tspan[2]-tspan[1])
function LindbladProp(H₀::Operator, Cₘ::Tuple{Vararg{Operator}}, Δt::Real)
    dimsmatch(H₀,Cₘ)
    F = real(eltype(H₀)); d = prod(dims(H₀))
    Id = Matrix{F}(I, d, d)
    # Constant Hamiltonian term
    L₀Δt = -1im*Δt.*(Id⊗full(H₀) .- transpose(full(H₀))⊗Id)
    # Add constant collapse operator terms
    L₀Δt .+= sum_collapse(Cₘ,Id,Δt)
    # Build constant propagator
    U = exp!(L₀Δt)
    return Propagator(U,Δt,dims(H₀))
end
LindbladProp(H₀::Operator, Hₙ::Tuple, Cₘ::Operator, tspan, steps::Integer) = LindbladProp(H₀,(Hₙ,),(Cₘ,),tspan,steps)
LindbladProp(H₀::Operator, Hₙ::Tuple, Cₘ::Tuple{Vararg{Operator}}, tspan, steps::Integer) = LindbladProp(H₀,(Hₙ,),Cₘ,tspan,steps)
LindbladProp(H₀::Operator, Hₙ::Tuple{Vararg{Tuple}}, Cₘ::Operator, tspan, steps::Integer) = LindbladProp(H₀,Hₙ,(Cₘ,),tspan,steps)
function LindbladProp(H₀::Operator, Hₙ::Tuple{Vararg{Tuple}}, Cₘ::Tuple{Vararg{Operator}}, tspan, steps::Integer)
    dimsmatch(H₀,first.(Hₙ)); dimsmatch(H₀,Cₘ)
    F = real(eltype(H₀)); d = prod(dims(H₀))
    Id = Matrix{F}(I, d, d)
    # Sampling times and spacing dt
    t₁, t₂ = tspan; dt = (t₂-t₁)/steps
    # Unpack constant and time dep operators
    H0, Hn = full(H₀), unpack_operators(1,full,Hₙ)
    # Build constant collapse propagator part
    U₀ = exp!(sum_collapse(Cₘ,Id,dt))
    # Multiply sampled propagators together to generate total evolution
    U = Matrix{Complex{F}}(I, size(U₀))
    H = Hermitian(zeros(compute_H_type(H0,Hn),size(H0)...))
    A = Matrix{Complex{F}}(undef,size(H0)...); B = similar(U); C = similar(H.data); D = similar(A)
    Λ = Vector{F}(undef,size(H,1))
    for i = 1:steps
        step_hamiltonian!(H.data,H0,Hn,(t₁,dt,i))
        expim!(A,H,Λ,C,D) # A = exp(-1im*H*dt)
        invA = Compat.LinearAlgebra.inv!(luf(A))
        mul!(B,U₀,U) # use the Lie product formula here for better results
        @static if VERSION < v"0.7.0-"
            At_mul_B!(U,invA⊗Id,B)
        else
            mul!(U,transpose(invA⊗Id),B)
        end
        I_kron_mul!(B,A,U)
        U,B = B,U # swappitty swap for the next step
    end
    return Propagator(U,float(t₂-t₁),dims(H₀))
end
LindbladProp(H) = throw(ArgumentError("invalid Propagator specification"))


# Utils
function unpack_operators(x::Number,f::Function,t::Tuple{Tuple{Operator,Function,Array},Vararg{Tuple}})
    a, b = first(t), unpack_operators(x,f,tail(t))
    return (x*f(a[1]), b[1]...), (a[2], b[2]...), (a[3], b[3]...)
end
function unpack_operators(x::Number,f::Function,t::Tuple{Tuple{Operator,Function},Vararg{Tuple}})
    a, b = first(t), unpack_operators(x,f,tail(t))
    return (x*f(a[1]), b[1]...), (a[2], b[2]...), ([], b[3]...)
end
unpack_operators(x::Number,f::Function,t::Tuple{}) = (), (), ()

function sum_collapse(Cₘ,Id,dt)
    # sum the Lindblad superoperators corresponding to the action of
    # D[Cᵢ](ρ) = CᵢρCᵢ† - 1/2 * (Cᵢ†Cᵢρ + ρCᵢ†Cᵢ)
    mapreduce(+,Cₘ) do Cᵢ
        C = full(Cᵢ); CdC = C'*C
        dt*(conj(C)⊗C - 0.5*(Id⊗CdC + transpose(CdC)⊗Id))
    end
end

function step_hamiltonian!(H,H0,Hn,tspec)
    Hₙ,fₙ,pₙ = Hn
    t₁, dt, i = tspec; tᵢ = t₁ + (i-1)*dt
    H .= -dt.*H0 # constant part
    for j = 1:length(Hₙ) # sample function at 3 points
        Hᵢ,fᵢ,pᵢ = Hₙ[j],fₙ[j],pₙ[j]
        H .+= -dt*(fᵢ(tᵢ,pᵢ) + fᵢ(tᵢ+0.5dt,pᵢ) + fᵢ(tᵢ+dt,pᵢ))/3 .* Hᵢ
    end
end

function compute_H_type(H0,Hn)
    T1 = promote_eltype(H0,Hn[1]...)
    T2 = promote_type(_typeof_H_funs(Hn[2],Hn[3])...)
    return promote_type(T1,T2)
end

@inline function _typeof_H_funs(funs::NTuple{N,Function},ps::NTuple{N,Vector}) where N
    this_T = typeof(first(funs)(0,first(ps)))
    return (this_T, _typeof_H_funs(tail(funs),tail(ps))...)
end

@inline _typeof_H_funs(funs::Tuple{},ps::Tuple{}) = ()
