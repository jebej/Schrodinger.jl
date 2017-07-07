"""
    expect(σ,ψ), expect(σ,ρ)

Compute the expectation value of an operator \$σ\$ given a state ket \$|ψ⟩\$ or a density matrix \$ρ\$. The expectation value is defined as
```math
\\begin{align*}
E(σ,|ψ⟩) &= ⟨ψ|σ|ψ⟩, \\\\
E(σ,ρ) &= \\textrm{tr}(σρ).
\\end{align*}
```
"""
expect(σ::Operator,ψ::Ket) = (dimsmatch(σ,ψ); dot(data(ψ),data(σ)*data(ψ)))
expect(ψ::Ket,σ::Operator) = expect(σ,ψ)
expect(σ::Operator,ρ::Operator) = (dimsmatch(σ,ρ); trace(data(σ)*data(ρ)))

# Faster expectation value mathods when evaluating many (see calc_expvals)
function expect!{T,D}(res,σ::Operator,states::Vector{Ket{T,D}})
    dimsmatch(σ,states[1])
    tmp = Vector{Complex128}(prod(dims(σ)))
    for (i,ψ) in enumerate(states)
        spmdv_mul!(tmp,data(σ),data(ψ))
        res[i] = dot(data(ψ),tmp)
    end
    return res
end

function expect!{T,D}(res,σ::Operator,states::Vector{Operator{T,D}})
    dimsmatch(σ,states[1])
    superσ = super(data(σ))
    tmp = Vector{Complex128}(prod(dims(σ))^2)
    N = length(tmp); sqrtNp1 = isqrt(N)+1
    fill!(res, zero(eltype(tmp)))
    for (i,ρ) in enumerate(states)
        spmdv_mul!(tmp,superσ,vec(data(ρ)))
        for n = 1:sqrtNp1:N
            res[i] += tmp[n]
        end
    end
    return res
end


levelprobs{T}(ψ::Ket{T,1}) = abs2.(data(ψ))
levelprobs{T,N}(ψ::Ket{T,N},s::Int) = diag(ptrace(ψ,ntuple_sans_m(s,Val{N})))
levelprobs{T,N}(ψ::Ket{T,N},out::Vector{Int}) = diag(ptrace(ψ,out))

function levelprobs{T,M}(states::Vector{Ket{T,M}},S=1:M)
    N = length(states)
    D = dims(states[1])
    probs = map(S) do s
        P = Matrix{Float64}(D[s],N)
        for n = 1:N
            P[:,n] = levelprobs(states[n],s)
        end
        return P.'
    end
    return probs
end

function levelprobs{T}(states::Vector{Ket{T,1}})
    N = length(states)
    P = Matrix{Float64}(dims(states[1])[1],N)
    for n = 1:N
        P[:,n] = levelprobs(states[n])
    end
    return P.'
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


function oper(vecA::AbstractVector,N=isqrt(length(vecA)))
    return reshape(vecA,(N,N))
end

function super(A::AbstractMatrix,B::AbstractMatrix=data(qeye(size(A,1))))
    return B.' ⊗ A
end
