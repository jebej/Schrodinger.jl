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

function gate_fidelity_choi(C::Operator,U::Operator)
    # calculate the average gate fidelity given a Choi matrix and a unitary gate
    length(dims(C)) == 2 && dims(C)[1]==dims(C)[2]==dims(U)[1] || throw(DimensionMismatch())
    V = qeye(dims(U)[1]) ⊗ U
    return gate_fidelity_choi(V'*C*V)
end

function gate_fidelity_choi(C::Operator)
    # calculate the average gate fidelity given a Choi matrix
    # Johnston, N. & Kribs, D. W. Quantum gate fidelity in terms of Choi
    # matrices. J. Phys. A Math. Theor. 44, 495303 (2011).
    D = dims(C)
    length(D)==2 && D[1]==D[2] || throw(ArgumentError("not a Choi matrix!"))
    d = D[1]
    return (d + real(sum(C[(i,i),(j,j)] for i=0:d-1,j=0:d-1)))/(d^2 + d)
end

function gate_fidelity_kraus(As::Vector{<:Operator})
    # calculate the average gate fidelity given a list of Kraus operators
    d = dims(As[1])[1]
    return (d + sum(abs2∘trace,As))/(d^2 + d)
end

function gate_fidelity_kraus(As::Vector{<:Operator},U::Operator)
    # calculate the average gate fidelity given a list of Kraus operators and
    # a unitary gate
    dimsmatch(As[1],U)
    d = dims(U)[1]
    return (d + sum(abs2∘inner,Base.product(As,(U,))))/(d^2 + d)
end

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
