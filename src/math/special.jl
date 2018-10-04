"""
    expect(σ,ψ), expect(σ,ρ)

Compute the expectation value of an operator \$σ\$ given a state ket \$|ψ⟩\$ or a density matrix \$ρ\$. The expectation value is defined as
```math
\\begin{align*}
E(σ,|ψ⟩) &= ⟨ψ|σ|ψ⟩, \\\\
E(σ,ρ) &= \\textrm{tr}(σρ).
\\end{align*}
```

A specialized method exists for vector of `Ket` or `Operator` inputs.
"""
expect(σ::Operator,ψ::Ket) = (dimsmatch(σ,ψ); dot(data(ψ),data(σ)*data(ψ)))
expect(ψ::Ket,σ::Operator) = expect(σ,ψ)
expect(σ::Operator,ρ::Operator) = (dimsmatch(σ,ρ); trace(data(σ)*data(ρ)))

# Faster expectation value methods when passing in many inputs
function expect!(res,σ::Operator,states::Vector{Ket{T,D}}) where {T,D}
    dimsmatch(σ,states[1])
    tmp = Vector{ComplexF64}(undef,prod(dims(σ)))
    for (i,ψ) in enumerate(states)
        mul!(tmp,data(σ),data(ψ))
        res[i] = dot(data(ψ),tmp)
    end
    return res
end

function expect!(res,σ::Operator,states::Vector{Operator{T,D}}) where {T,D}
    dimsmatch(σ,states[1])
    superσ = super(data(σ))
    tmp = Vector{ComplexF64}(undef,prod(dims(σ))^2)
    N = length(tmp); sqrtNp1 = isqrt(N)+1
    fill!(res, zero(eltype(tmp)))
    for (i,ρ) in enumerate(states)
        mul!(tmp,superσ,vec(data(ρ)))
        for n = 1:sqrtNp1:N
            res[i] += tmp[n]
        end
    end
    return res
end

"""
    levelprobs(ψ), levelprobs(ψ,s)

Compute the level occupation probabilities. For a `Ket`, this simply corresponds to the absolute square of the amplitude of each level. For an `Operator`, the function returns the diagonal.

A system index, or vector of indices, can be passed as a second argument. In that case, the full system will first be partial traced to keep only the desired index. Level occupation probabilities are then calculated from the resulting reduced density matrix. If a vector of indices is passed, occupation probabilities are calculated for a fully reduced density matrix for each index.

A specialized method exists for vector of `Ket` or `Operator` inputs.
"""
levelprobs(ψ::Ket{T,N}) where {T,N} = abs2(ψ)
levelprobs(ψ::Ket{T,N},s::Integer) where {T,N} = real(diag(ptrace(ψ,ntuple_sans_m(s,Val(N)))))
levelprobs(ρ::Operator{T,N}) where {T,N} = real(diag(ρ))
levelprobs(ρ::Operator{T,N},s::Integer) where {T,N} = real(diag(ptrace(ρ,ntuple_sans_m(s,Val(N)))))
levelprobs(ψ::QuObject,S::AbstractVector) = map(s->levelprobs(ψ,s),S)

# Faster levelprobs methods when passing in many inputs
function levelprobs(states::Vector{Ket{T,M}},S::Union{Integer,AbstractVector}) where {T,M}
    N = length(states)
    D = dims(states[1])
    probs = map(S) do s
        s > M && throw(ArgumentError("System index $s is larger than the number of systems, $M."))
        P = Matrix{Float64}(undef,D[s],N)
        for n = 1:N
            P[:,n] = levelprobs(states[n],s)
        end
        return transpose(P)
    end
    return probs
end

function levelprobs(states::Vector{Ket{T,M}}) where {T,M}
    N = length(states)
    P = Matrix{Float64}(undef,prod(dims(states[1])),N)
    for n = 1:N
        P[:,n] = levelprobs(states[n])
    end
    return transpose(P)
end

"""
    fidelity(ρ,σ), fidelity(ρ,ψ), fidelity(ψ,ϕ)

Compute the fidelity between density matrices \$ρ\$ and \$σ\$, a density matrix \$ρ\$ and a ket \$|ψ⟩\$, or two kets \$|ψ⟩\$ and \$|ϕ⟩\$. The fidelity in those three cases is defined as
```math
\\begin{align*}
F(ρ,σ) &= \\textrm{tr}\\sqrt{ρ^{1/2}σρ^{1/2}}, \\\\
F(ρ,|ψ⟩) &= \\sqrt{⟨ψ|ρ|ψ⟩}, \\\\
F(|ψ⟩,|ϕ⟩) &= \\left|⟨ψ|ϕ⟩\\right|.
\\end{align*}
```
See also [`fidelity2`](@ref), which is the square of the state fidelity.
"""
function fidelity(ρ::Operator,σ::Operator)
    dimsmatch(ρ,σ)
    ishermitian(ρ)&&ishermitian(σ) || throw(ArgumentError("the operators must be Hermitian"))
    sqrtρ = sqrtm(Hermitian(full(ρ)))
    A = sqrtρ*data(σ)*sqrtρ
    D = eigvals(A)
    res = 0.0
    for i = 1:prod(dims(ρ))
        d = real(D[i])
        (d>0) && (res += sqrt(d))
    end
    return res
end
fidelity(ρ::Operator,ψ::Ket) = sqrt(fidelity2(ρ,ψ))
fidelity(ψ::Ket,ρ::Operator) = fidelity(ρ,ψ)
fidelity(ψ::Ket,ϕ::Ket) = abs(dot(ψ,ϕ))

"""
    fidelity2(ρ,ψ)

Compute the Uhlmann state fidelity between density matrices \$ρ\$ and \$σ\$, a density matrix \$ρ\$ and a ket \$|ψ⟩\$, or two kets \$|ψ⟩\$ and \$|ϕ⟩\$. The Uhlmann state fidelity is simply defined as the square of the "regular" state [`fidelity`](@ref).
"""
fidelity2(ρ::Operator,σ::Operator) = abs2(fidelity(ρ,σ))
function fidelity2(ρ::Operator,ψ::Ket)
    ishermitian(ρ) || throw(ArgumentError("the operator must be Hermitian"))
    return abs(expect(ρ,ψ))
end
fidelity2(ψ::Ket,ρ::Operator) = fidelity2(ρ,ψ)
fidelity2(ψ::Ket,ϕ::Ket) = abs2(dot(ψ,ϕ))

function fidelity2(A::NTuple{N,T},ϕ::Ket) where {N,T<:Ket}
    # If sum(A) = ψ is a valid Ket, then this function will calculate the state fidelity between ψ and ϕ, ignoring all relative phases between the different parts of ψ (as well as the global phase, of course)
    # F² = ∑ᵢ(|⟨ψᵢ,ϕ⟩|² + ∑ⱼ₌₁ⁱ⁻¹2*|⟨ψᵢ,ϕ⟩|*|⟨ψⱼ,ϕ⟩|)
    f2 = fidelity2.(A,ϕ) # (|⟨ψ₁,ϕ⟩|², |⟨ψ₂,ϕ⟩|², ...)
    res = 0.0
    for i = 1:length(f2)
        res += f2[i]
        for j = 1:i-1
            res += 2*sqrt(f2[i])*sqrt(f2[j])
        end
    end
    return res
end


"""
    entanglement_fidelity(U,V)

Compute the entanglement fidelity of an imperfect quantum operation `V` with respect to an ideal operation `U`. The entanglement fidelity is given by the following formula:
```math
F_\\text{e}(\\mathcal{E}) = ⟨ϕ|(I⊗\\mathcal{E})(ϕ)|ϕ⟩,
```
where \$|ϕ⟩ = \\frac{1}{\\sqrt{d}}\\sum_{i=1}^d|i⟩|i⟩\$ is the maximally entangled state and \$\\mathcal{E}\$ is a trace preserving quantum operation.

When trying to find the fidelity of a quantum gate \$U\$, implemented imperfectly as \$V\$, then \$\\mathcal{E} = U^†∘V\$. In the case where \$U\$ and \$V\$ are a simple operators,
```math
\\begin{align*}
(I⊗\\mathcal{E})(ϕ) = (I⊗U^†V)|ϕ⟩⟨ϕ|(I⊗U^†V)^† \\quad\\text{and so}
\\end{align*}
```
```math
\\begin{align*}
F_\\text{e}(\\mathcal{E}) &= ⟨ϕ|(I⊗U^†V)|ϕ⟩⟨ϕ|(I⊗U^†V)^†|ϕ⟩ \\\\
  &= |⟨ϕ|(I⊗U^†V)|ϕ⟩|^2 \\\\
  &= |⟨U,V⟩|^2/d^2 \\\\
\\end{align*}
```

The entanglement fidelity is related to the [`average gate fidelity`](@ref gate_fidelity) by
```math
\\begin{align*}
\\overline{F}(\\mathcal{E}) =  \\frac{d F_\\text{e}(\\mathcal{E}) + 1}{d + 1}.
\\end{align*}
```
"""
function entanglement_fidelity(U::Operator,V::Operator,d::Integer=prod(dims(U)))
    dimsmatch(U,V)
    return abs2(inner(U,V))/d^2
end

"""
    gate_fidelity(U,V)

Compute the average gate fidelity of an imperfect quantum operation `V` with respect to an ideal operation `U`.

The average gate fidelity is related to the [`entanglement fidelity`](@ref entanglement_fidelity) by
```math
\\begin{align*}
\\overline{F}(\\mathcal{E}) =  \\frac{d F_\\text{e}(\\mathcal{E}) + 1}{d + 1}.
\\end{align*}
```
"""
function gate_fidelity(U::Operator,V::Operator,d::Integer=prod(dims(U)))
    dimsmatch(U,V)
    return (abs2(inner(U,V)) + d)/(d^2 + d)
end

function oper(vecA::AbstractVector,N=isqrt(length(vecA)))
    return reshape(vecA,(N,N))
end

function super(A::AbstractMatrix,B::AbstractMatrix=data(qeye(size(A,1))))
    return transpose(B) ⊗ A
end
