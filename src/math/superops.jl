# assumes Choi matrix
apply_process(C::Operator,ψ::Ket) = apply_process(C,Operator(ψ))

apply_process(C::Operator,ρ::Operator{T,D}) where {T,D} =
    ptrace((transpose(ρ)⊗qeye(dims(ρ)))*C,ntuple(identity,D))

# assumes Kraus operators
apply_process(As::Vector{<:Operator},ρ::Operator) = sum(A->A*ρ*A',As)

function gate_fidelity_choi(C::Operator,U::Operator)
    # calculate the average gate fidelity given a Choi matrix and a unitary gate
    dims(C) != (dims(U)...,dims(U)...) && throw(DimensionMismatch())
    V = qeye(dims(U)) ⊗ U
    return gate_fidelity_choi(V'*C*V)
end

function gate_fidelity_choi(C::Operator)
    # calculate the average gate fidelity given a Choi matrix
    # Johnston, N. & Kribs, D. W. Quantum gate fidelity in terms of Choi
    # matrices. J. Phys. A Math. Theor. 44, 495303 (2011).
    d = check_choi(C)
    TrA² = let CC = Operator(data(C),(d,d)) # rewrap operator to use tensored indexing
        real(sum(CC[(i,i),(j,j)] for i=0:d-1, j=0:d-1))
    end
    return (d + TrA²)/(d^2 + d)
end

function gate_fidelity_kraus(As::Vector{Operator{T,D}}) where {T,D}
    # calculate the average gate fidelity given a list of Kraus operators
    d = prod(i->dims(first(As))[i],1:D)
    return (d + sum(abs2∘trace,As))/(d^2 + d)
end

function gate_fidelity_kraus(As::Vector{<:Operator{T1,D}},U::Operator{T2,D}) where {T1,T2,D}
    # calculate the average gate fidelity given a list of Kraus operators and
    # a unitary gate
    dimsmatch(As,U)
    d = prod(i->dims(U)[i],1:D)
    return (d + sum(abs2∘inner,product(As,(U,))))/(d^2 + d)
end

# conversion functions
"""
    operator_to_choi(O)

Compute the Choi matrix ``C_O`` of a unitary map represented by an operator `O`. We use the
convention of applying `O` to the second half of the extended Hilbert space, i.e.:
```math
C_O = ∑_{i,j} |i⟩⟨j| ⊗ O|i⟩⟨j|O^†
```
"""
operator_to_choi(O::Operator) = Operator(vec(O))

kraus_to_choi(As::Vector{<:Operator}) = sum(operator_to_choi,As)

kraus_to_natural(As::Vector{<:Operator}) = sum(A->conj(A)⊗A,As)

# go through Choi representation; the extra step is cheap compared to the eigendecomposition
natural_to_kraus(K::Operator) = choi_to_kraus(natural_to_choi(K))

# the conversion between the natural and Choi is a self-inverse
natural_to_choi(K::Operator) = choi_to_natural(K)

function choi_to_natural(C::Operator)
    # the "natural" representation is the same as Choi with shuffled dimensions
    d = check_choi(C)
    K = reshape(permutedims(reshape(data(C),(d,d,d,d)),[1,3,2,4]),(d^2,d^2))
    return Operator(K,dims(C))
end

function choi_to_kraus(C::Operator)
    check_choi(C)
    D,V = eigen(Hermitian(data(C)),1E-10,Inf)
    dm = _half_dims(dims(C))
    return [Operator(unvec(√(D[i])*V[:,i]),dm) for i = length(D):-1:1]
end

function choi_to_chi(C::Operator{T,D}) where {T,D}
    B = Operator(_choi_to_pauli_basis(_half_dims(dims(C))),dims(C))
    return (B*C*B')/2^(D÷2)
end

function chi_to_choi(χ::Operator{T,D}) where {T,D}
    B = Operator(_choi_to_pauli_basis(_half_dims(dims(χ))),dims(χ))
    return (B'*χ*B)/2^(D÷2)
end

function _choi_to_pauli_basis(D::Dims{N}) where N
    all(d->d==2,D) || throw(ArgumentError("only valid with qubits!"))
    B = Matrix{Float64}(undef,4^N,4^N)
    # change to I,X,Y,Z basis (Y = -i*σy)
    Ps = ([1 0; 0 1], [0 1; 1 0], [0 -1; 1 0], [1 0; 0 -1])
    for (i,ops) in enumerate(product(ntuple(_->Ps,N)...))
        B[:,i] = vec(reduce(⊗,ops))
    end
    return B
end

Base.vec(O::Operator) = Ket(vec(data(O)),(dims(O)...,dims(O)...))

unvec(v::Ket) = Operator(unvec(data(v)),_half_dims(dims(v)))

super(A::Operator) = conj(A) ⊗ A # represents A*ρ*A†
super(A::Operator,B::Operator) = transpose(B) ⊗ A

# linear operator that represents the action AρB on a vectorized version of ρ
super(A::AbstractMatrix,B::AbstractMatrix=data(qeye(size(A,1)))) = transpose(B) ⊗ A

function unvec(vecA::AbstractVector)
    # unvectorize a vector into a square matrix
    d = isqrt(length(vecA))
    return reshape(vecA,(d,d))
end

function check_choi(C::Operator{T,D}) where {T,D}
    dm = dims(C)
    if D%2 != 0 || prod(i->dm[i],1:D÷2) != prod(i->dm[i],D÷2+1:D)
        throw(ArgumentError("not a valid Choi matrix!"))
    end
    return prod(i->dm[i],1:D÷2)
end

_half_dims(dm::Dims{D}) where{D} = ntuple(i->dm[i],D÷2)
