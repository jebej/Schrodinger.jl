# Schrodinger equation
SchrodingerEvo(H₀::Operator) = Liouvillian(dims(H₀),-1im*data(H₀))

function SchrodingerEvo(H₀::Operator, Hₙ::Tuple{Vararg{Tuple}})
    dimsmatch(H₀,first.(Hₙ))
    # Unpack time-dependent Hamiltonian terms
    Lₙ,fₙ,pₙ = unpack_operators(Hₙ)
    return Liouvillian(dims(H₀),-1im*data(H₀),-1im.*data.(Lₙ),fₙ,pₙ)
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
    L₀ = -1im*super_vonneumann(H₀,Id)
    # Add constant collapse operator terms
    L₀ += sum_collapse(Cₘ,Id)
    return Liouvillian(dims(H₀),L₀)
end
function LindbladEvo(H₀::Operator, Hₙ::Tuple{Vararg{Tuple}}, Cₘ::Tuple{Vararg{Operator}})
    dimsmatch(H₀,first.(Hₙ)); dimsmatch(H₀,Cₘ)
    Id = data(qeye(dims(H₀)))
    # Constant Hamiltonian term
    L₀ = -1im*super_vonneumann(H₀,Id)
    # Add constant collapse operator terms
    L₀ += sum_collapse(Cₘ,Id)
    # Unpack time-dependent Hamiltonian terms
    Lₙ,fₙ,pₙ = unpack_operators(Hₙ)
    return Liouvillian(dims(H₀),L₀,-1im.*super_vonneumann.(Lₙ,(Id,)),fₙ,pₙ)
end
function LindbladEvo(H₀::Operator, Hₙ::Tuple{Vararg{Tuple}}, Cₘ::Tuple{Vararg{Tuple}})
    dimsmatch(H₀,first.(Hₙ)); dimsmatch(H₀,first.(Cₘ))
    Id = data(qeye(dims(H₀)))
    # Constant Hamiltonian term
    L₀ = -1im*super_vonneumann(H₀,Id)
    # Unpack time-dependent Hamiltonian terms
    Lₙ,fₙ,pₙ = unpack_operators(Hₙ)
    # Unpack time-dependent collapse operator terms
    Lₘ,fₘ,pₘ = unpack_operators(Cₘ)
    # Superoperatorize time-dependent operators
    Lₙ_super = -1im.*super_vonneumann.(Lₙ,(Id,))
    Lₘ_super = super_collapse.(data.(Lₘ),(Id,))
    # Merge tuples
    return Liouvillian(dims(H₀),L₀,(Lₙ_super...,Lₘ_super...),(fₙ...,fₘ...),(pₙ...,pₘ...))
end
LindbladEvo(H₀::Operator,Cₘ::Vararg{Operator}) = LindbladEvo(H₀,Cₘ)
LindbladEvo(H::Tuple{Operator,Vararg{Tuple}},Cₘ::Vararg{Operator}) = LindbladEvo(first(H),tail(H),Cₘ)
LindbladEvo(H::Tuple{Operator,Vararg{Tuple}},Cₘ::Tuple{Vararg{Operator}}) = LindbladEvo(first(H),tail(H),Cₘ)
LindbladEvo(::Any) = throw(ArgumentError("invalid Hamiltonian specification"))


# Schrodinger evo propagator
SchrodingerProp(H₀::Operator, tspan) = SchrodingerProp(H₀,float(tspan[2]-tspan[1]))
function SchrodingerProp(H₀::Operator, Δt::Real)
    # Constant Hamiltonian term
    U = expim(Hermitian(-Array(H₀).*Δt))
    return Propagator(U,Δt,dims(H₀))
end
SchrodingerProp(H₀::Operator, Hₙ::Tuple, tspan, steps::Integer) = SchrodingerProp(H₀,(Hₙ,),tspan,steps)
function SchrodingerProp(H₀::Operator, Hₙ::Tuple{Vararg{Tuple}}, tspan, steps::Integer)
    dimsmatch(H₀,first.(Hₙ)); F = real(eltype(H₀))
    # Sampling times and spacing dt
    t₁, t₂ = tspan; dt = (t₂-t₁)/steps
    # Unpack time-dependent Hamiltonian terms
    H_tdep = unpack_operators(Hₙ)
    # Multiply sampled propagators together to generate total evolution
    U = Matrix{Complex{F}}(I, size(H₀))
    H = Hermitian(zeros(compute_H_type(H₀,H_tdep),size(H₀)...))
    A = similar(U); B = similar(U); C = similar(U)
    for i = 1:steps
        step_hamiltonian!(H.data,H₀,H_tdep,(t₁,dt,i)) # calc H for this time step
        expim!(A,H,C) # A = exp(-1im*H*dt)
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
    Id = Diagonal(I, prod(dims(H₀)))
    # Constant Hamiltonian term
    L₀Δt = (-1im*Δt).*(Id⊗Array(H₀) .- transpose(Array(H₀))⊗Id)
    # Add constant collapse operator terms
    L₀Δt .+= Δt.*sum_collapse(dense.(Cₘ),Id)
    # Build constant propagator
    U = LinearAlgebra.exp!(L₀Δt)
    return Propagator(U,Δt,dims(H₀))
end
LindbladProp(H₀::Operator, Hₙ::Tuple, Cₘ::Operator, tspan, steps::Integer) = LindbladProp(H₀,(Hₙ,),(Cₘ,),tspan,steps)
LindbladProp(H₀::Operator, Hₙ::Tuple, Cₘ::Tuple{Vararg{Operator}}, tspan, steps::Integer) = LindbladProp(H₀,(Hₙ,),Cₘ,tspan,steps)
LindbladProp(H₀::Operator, Hₙ::Tuple{Vararg{Tuple}}, Cₘ::Operator, tspan, steps::Integer) = LindbladProp(H₀,Hₙ,(Cₘ,),tspan,steps)
function LindbladProp(H₀::Operator, Hₙ::Tuple{Vararg{Tuple}}, Cₘ::Tuple{Vararg{Operator}}, tspan, steps::Integer)
    dimsmatch(H₀,first.(Hₙ)); dimsmatch(H₀,Cₘ); F = real(eltype(H₀));
    Id = Diagonal(I, prod(dims(H₀)))
    # Sampling times and spacing dt
    t₁, t₂ = tspan; dt = (t₂-t₁)/steps
    # Unpack time-dependent Hamiltonian terms
    H_tdep = unpack_operators(Hₙ)
    # Build constant collapse propagator part
    U₀ = LinearAlgebra.exp!(dt.*sum_collapse(dense.(Cₘ),Id))
    # Multiply sampled propagators together to generate total evolution
    U = Matrix{Complex{F}}(I, size(U₀))
    H = Hermitian(zeros(compute_H_type(H₀,H_tdep),size(H₀)...))
    A = Matrix{Complex{F}}(undef,size(H₀)...); B = similar(U); C = similar(A)
    for i = 1:steps
        step_hamiltonian!(H.data,H₀,H_tdep,(t₁,dt,i))
        expim!(A,H,C) # A = exp(-1im*H*dt)
        invA = LinearAlgebra.inv!(lu(A))
        #
        mul!(B,U₀,U) # TODO: use the Lie product formula here for better results
        mul!(U,transpose(invA⊗Id),B)
        I_kron_mul!(B,A,U)
        U,B = B,U # swappitty swap for the next step
    end
    return Propagator(U,float(t₂-t₁),dims(H₀))
end
LindbladProp(H) = throw(ArgumentError("invalid Propagator specification"))


# Utils
function unpack_operators(t::Tuple{Tuple{Operator,Function,Array},Vararg{Tuple}})
    a, b = first(t), unpack_operators(tail(t))
    return (a[1], b[1]...), (a[2], b[2]...), (a[3], b[3]...)
end
function unpack_operators(t::Tuple{Tuple{Operator,Function},Vararg{Tuple}})
    a, b = first(t), unpack_operators(tail(t))
    return (a[1], b[1]...), (a[2], b[2]...), (Float64[], b[3]...)
end
unpack_operators(t::Tuple{}) = (), (), ()

function super_vonneumann(H,Id)
    # represents the action [H,ρ] as a linear superoperator I⊗H - Hᵀ⊗I
    return Id⊗data(H) - copy(transpose(data(H)))⊗Id
end

function sum_collapse(Cₘ, Id)
    # sum the Lindblad collapse superoperators
    mapreduce(Cᵢ->super_collapse(data(Cᵢ), Id), +, Cₘ)
end

function super_collapse(C::AbstractMatrix, Id::AbstractMatrix)
    # represent the Lindblad superoperator corresponding to the action of
    # D[C](ρ) = CρC† - 1/2 * (C†Cρ + ρC†C) as a linear superoperator
    CdC = C'*C
    return conj(C)⊗C - (Id⊗CdC + copy(transpose(CdC))⊗Id)./2
end

function step_hamiltonian!(H,H0,Hn,tspec)
    Hₙ,fₙ,pₙ = Hn
    t₁, dt, i = tspec; tᵢ = t₁ + (i-1)*dt
    H .= -dt.*data(H0) # constant part
    for j = 1:length(Hₙ) # sample function at 3 points
        Hᵢ,fᵢ,pᵢ = data(Hₙ[j]),fₙ[j],pₙ[j]
        H .+= (-dt*(fᵢ(tᵢ,pᵢ) + fᵢ(tᵢ+0.5dt,pᵢ) + fᵢ(tᵢ+dt,pᵢ))/3) .* Hᵢ
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
