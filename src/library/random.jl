"""
    rand_unitary(N, dims=(N,))

Generate a Haar distributed random unitary operator for a Hilbert space of size `N`. It is possible to specify the subspace dimensions with the `dims` argument. Returns a dense matrix.

# Example
```jldoctest
julia> U = rand_unitary(4,(2,2));

julia> U'*U ≈ qeye(4,(2,2))
true
```
"""
function rand_unitary(N::Integer, dims::Dims=(N,))
    return Operator(rand_unitary(ComplexF64,N),dims)
end

function rand_unitary(::Type{T},n::Integer) where {T<:LinearAlgebra.BlasComplex}
    # Generate a Haar distributed random unitary matrix
    # ref https://arxiv.org/pdf/math-ph/0609050.pdf
    A = randn(T,n,n)
    f = qr!(A)
    U = Matrix(LinearAlgebra.QRCompactWYQ(f.factors,f.T))
    @inbounds for i = 1:n
        U[:,i] .*= sign(f.factors[i,i]) # U = Q*Diagonal(diag(R)./abs.(diag(R)))
    end
    return U
end
