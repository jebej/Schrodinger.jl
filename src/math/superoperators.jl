apply_process(C::Operator,ψ::Ket) = apply_process(C,Operator(ψ))
apply_process(C::Operator,ρ::Operator) = ptrace((transpose(ρ)⊗qeye(size(ρ,1)))*C,1)

function operator_to_choi(O::Operator)
    length(dims(O)) == 1 || throw(ArgumentError("multi-space operators not supported yet!"))
    d = dims(O)[1]
    C = dense(qzero(eltype(O),d^2,(d,d)))
    for i=1:d, j=1:d, i′=1:d, j′=1:d
        @inbounds C[(i-1,i′-1),(j-1,j′-1)] = O[i,i′]*O[j,j′]'
    end
    return C
end

gate_fidelity(C) = (d=dims(C)[1]; (sum(@inbounds C[(i,i),(j,j)] for i=0:1,j=0:1)/d+1)/(d+1))

function choi_to_kraus(C)
    D,V = eigen(Hermitian(data(C)),1E-8,Inf)
    return [Operator(unvec(√(D[i])*V[:,i])) for i = length(D):-1:1]
end

function unvec(vecA::AbstractVector)
    # unvectorize a vector into a square matrix
    d = isqrt(length(vecA))
    return reshape(vecA,(d,d))
end

function super(A::AbstractMatrix,B::AbstractMatrix=data(qeye(size(A,1))))
    return transpose(B) ⊗ A
end
