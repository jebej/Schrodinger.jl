function expect(σ::Operator,ψ::Ket)
    dimsmatch(σ,ψ)
    return dot(data(ψ),data(σ)*data(ψ))
end
expect(ψ::Ket,σ::Operator) = expect(σ,ψ)

function expect(σ::Operator,ρ::Operator)
    dimsmatch(σ,ρ)
    return trace(data(σ)*data(ρ))
end

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


tuple_sans_m{n}(m,::Type{Val{n}}) = sorted_setdiff(ntuple(identity,Val{n}),(m,))

calc_probs{T}(ψ::Ket{T,1}) = abs2.(data(ψ))
calc_probs{T,N}(ψ::Ket{T,N},s::Int) = diag(ptrace(ψ,tuple_sans_m(s,Val{N})))
calc_probs{T,N}(ψ::Ket{T,N},out::Vector{Int}) = diag(ptrace(ψ,out))

function calc_probs{T,M}(states::Vector{Ket{T,M}})
    N = length(states)
    S = length(dims(states[1]))
    probs = map(1:S) do s
        d = dims(states[1])[s]
        P = Matrix{Float64}(d,N)
        out = [1:s-1;s+1:N]
        for n = 1:N
            P[:,n] = calc_probs(states[n],s)
        end
        return P.'
    end
    return probs
end

function calc_probs{T}(states::Vector{Ket{T,1}})
    N = length(states)
    d = dims(states[1])[1]
    P = Matrix{Float64}(d,N)
    for n = 1:N
        P[:,n] = calc_probs(states[n])
    end
    return P.'
end

function oper(vecA::AbstractVector,N=isqrt(length(vecA)))
    return reshape(vecA,(N,N))
end

function super(A::AbstractMatrix,B::AbstractMatrix=data(qeye(size(A,1))))
    return B.' ⊗ A
end
