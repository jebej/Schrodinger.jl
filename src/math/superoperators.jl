apply_process(C::Operator,ψ::Ket) = apply_process(C,Operator(ψ))
# assumes Choi matrix
apply_process(C::Operator,ρ::Operator{T,D}) where {T,D} =
    ptrace((transpose(ρ)⊗qeye(dims(ρ)))*C,ntuple(identity,D))
# assumes Kraus operators
apply_process(As::Vector{<:Operator},ρ::Operator) = sum(A->A*ρ*A',As)

function gate_fidelity_choi(C::Operator,U::Operator)
    # calculate the average gate fidelity given a Choi matrix and a unitary gate
    if dims(C) != (dims(U)...,dims(U)...)
        throw(DimensionMismatch())
    end
    V = qeye(dims(U)) ⊗ U
    return gate_fidelity_choi(V'*C*V)
end

function gate_fidelity_choi(C::Operator{T,D}) where {T,D}
    # calculate the average gate fidelity given a Choi matrix
    # Johnston, N. & Kribs, D. W. Quantum gate fidelity in terms of Choi
    # matrices. J. Phys. A Math. Theor. 44, 495303 (2011).
    dm = dims(C)
    if D%2 != 0 || prod(i->dm[i],1:D÷2) != prod(i->dm[i],D÷2+1:D)
        throw(ArgumentError("not a valid Choi matrix!"))
    end
    d = prod(i->dm[i],1:D÷2)
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

Compute the Choi-Jamilkowski matrix representation of an operator `O`. This uses
the convention of applying the action of `O` to the second half of the extended
Hilbert space, i.e.:

$$
J(\Phi) = (I ⊗ \Phi) |Ω><Ω|
$$
"""
function operator_to_choi(O::Operator)
    newdims = (dims(O)...,dims(O)...)
    return Operator(vec(data(O))*vec(data(O))',newdims)
end

function kraus_to_natural(As::Vector{<:Operator})
    return sum(A->conj(A)⊗A,As)
end

function natural_to_kraus(K::Operator)
    # we go through the Choi representation since the extra conversion step is
    # basically free compared to the eigendecomposition
    return choi_to_kraus(natural_to_choi(K))
end

function natural_to_choi(K::Operator)
    # the conversion between the natural and Choi is a self-inverse
    return choi_to_natural(K)
end

function choi_to_natural(C::Operator)
    # the "natural" representation is the same as Choi with shuffled dimensions
    length(dims(C)) == 2 || throw(ArgumentError("multi-space operators not supported yet!"))
    d = dims(C)[1]
    K = reshape(permutedims(reshape(data(C),(d,d,d,d)),[1,3,2,4]),(d^2,d^2))
    return Operator(K,dims(C))
end

function choi_to_kraus(C::Operator)
    length(dims(C)) == 2 || throw(ArgumentError("multi-space operators not supported yet!"))
    @static if VERSION < v"0.7.0-"
        @inbounds for i = 1:size(C,1); C[i,i] = real(C[i,i]); end
    end
    D,V = eigen(Hermitian(data(C)),1E-10,Inf)
    return [Operator(unvec(√(D[i])*V[:,i]),(dims(C)[1],)) for i = length(D):-1:1]
end

function kraus_to_choi(As::Vector{<:Operator})
    length(dims(As[1])) == 1 || throw(ArgumentError("multi-space operators not supported yet!"))
    d = dims(As[1])[1]
    return Operator(sum(A->vec(data(A))*vec(data(A))',As),(d,d))
end

function choi_to_chi(C::Operator)
    length(dims(C)) == 2 || throw(ArgumentError("multi-space operators not supported yet!"))
    B = Operator(_choi_to_pauli_basis((dims(C)[1],)),dims(C))
    return (B*C*B')/2^1 # TODO: fix for multi-qubit
end

function chi_to_choi(χ::Operator)
    length(dims(χ)) == 2 || throw(ArgumentError("multi-space operators not supported yet!"))
    B = Operator(_choi_to_pauli_basis((dims(χ)[1],)),dims(χ))
    return (B'*χ*B)/2^1 # TODO: fix for multi-qubit
end

function _choi_to_pauli_basis(D::NTuple{N,Int}) where N
    all(d->d==2,D) || throw(ArgumentError("only valid with qubits!"))
    B = Matrix{Float64}(undef,4^N,4^N)
    # change to I,X,Y,Z basis (Y = -i*σy)
    Ps = ([1 0; 0 1], [0 1; 1 0], [0 -1; 1 0], [1 0; 0 -1])
    for (i,ops) in enumerate(product(ntuple(_->Ps,N)...))
        B[:,i] = vec(reduce(⊗,ops))
    end
    return B
end

Base.vec(O::Operator) = Ket(copy(vec(data(O))),dims(O).^2)

unvec(v::Ket) = Operator(copy(unvec(data(O))),isqrt.(dims(v)))

super(A::Operator) = conj(A) ⊗ A # represents A*ρ*A†
super(A::Operator,B::Operator) = transpose(B) ⊗ A

function unvec(vecA::AbstractVector)
    # unvectorize a vector into a square matrix
    d = isqrt(length(vecA))
    return reshape(vecA,(d,d))
end

function super(A::AbstractMatrix,B::AbstractMatrix=data(qeye(size(A,1))))
    # represents the action AρB, on a vectorized version of ρ
    return transpose(B) ⊗ A
end
