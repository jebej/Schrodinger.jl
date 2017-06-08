"""
    basis(N, n)

Generate a basis state (a.k.a. Fock or number state) ket \$|n⟩\$, in a Hilbert space of size `N`. Note that the size of the Hilbert space must be *at least* `n+1`. The function `fock` is an alias for `basis`.

Returns a sparse vector.

# Example
```jldoctest
julia> ψ = basis(3,2)
3-d Schrodinger.Ket{SparseVector{Float64,Int64},1} with space dimensions 3:
1.00∠0°|2⟩
```
"""
function basis(N::Integer, n::Integer)
    N>n || throw(ArgumentError("basis level $n is too large for a $N-d space"))
    # Julia is 1-indexed!
    return Ket(SparseVector(N,[n+1],[1.0]),(N,))
end


"""
    coherent(N, α, analytic=false)

Generate a coherent state ket \$|α⟩\$, in a Hilbert space of size `N`. To create a coherent density operator, use the `Operator` function: `Operator(coherent(N,n))`.

Two methods can be used for generating a coherent state: via application of a displacment operator on a ground state (the default), or analytically, with the formula

```math
|α⟩ = e^{-\\frac{|α|^2}{2}} \\sum_{n=0}^{N-1} \\frac{α^n}{\\sqrt{n!}} |n⟩.
```

While the operator method will return a normalized ket, the analytic method will not. Both methods converge as `N` gets larger. The analytic method is also much faster, especially for large `N`.

Returns a dense vector.

# Example
```jldoctest
julia> coherent(6,0.4+1im)
6-d Schrodinger.Ket{Array{Complex{Float64},1},1} with space dimensions 6:
0.60∠68°|1⟩ + 0.56∠0°|0⟩ + 0.46∠136°|2⟩ + 0.29∠-155°|3⟩ + 0.15∠-87°|4⟩
```
"""
function coherent(N::Integer, α::Number, analytic::Bool=false)
    if analytic let a = exp(-abs2(α)/2) # due to julia issue #15276
        x = [a*α^n/sqrtfact(n) for n = 0:N-1]
        return Ket(x,(N,))
    end else
        return Ket(data(displacementop(N,α))[:,1],(N,)) # first column of D(α)
    end
end

"""
    thermal(N, n)

Generate a thermal state density matrix \$ρ_n\$ with particle number `n`, in a Hilbert space of size `N`. A thermal state \$ρ_n\$ is a probabilistic mixture of basis states such that the expectation value of the number operator \$\\hat{n}\$ is `n`. Note that this is true only if \$N≫n\$. The returned density matrix is always normalized.

Returns a sparse matrix.

# Example
```jldoctest
julia> N=5; n=0.2;

julia> ρ = thermal(N,n)
5×5 Schrodinger.Operator{Schrodinger.Herm,SparseMatrixCSC{Float64,Int64},1} with space dimensions 5:
 0.833441  0.0       0.0        0.0         0.0
 0.0       0.138907  0.0        0.0         0.0
 0.0       0.0       0.0231511  0.0         0.0
 0.0       0.0       0.0        0.00385852  0.0
 0.0       0.0       0.0        0.0         0.000643087

julia> expect(numberop(N),ρ)
0.19935691318327978
```
"""
function thermal(N::Integer, n::Real)
    β = log(1.0/n+1.0)
    rowval = collect(1:N)
    colptr = Vector{Int}(N+1); colptr[1:N] = rowval; colptr[end] = N+1
    nzval  = normalize!([exp(-β*k) for k = 0:N-1],1)
    return Operator(SparseMatrixCSC(N,N,colptr,rowval,nzval),(N,))
end

"""
    maxmixed(N)

Generate a maximally mixed density matrix. The maximally mixed state is a mixture of basis states with uniform probability.

Returns a sparse matrix.

# Example
```jldoctest
julia> maxmixed(4)
4×4 Schrodinger.Operator{Schrodinger.Herm,SparseMatrixCSC{Float64,Int64},1} with space dimensions 4:
 0.25  0.0   0.0   0.0
 0.0   0.25  0.0   0.0
 0.0   0.0   0.25  0.0
 0.0   0.0   0.0   0.25
```
"""
function maxmixed(N::Integer)
    rowval = collect(1:N)
    colptr = Vector{Int}(N+1); colptr[1:N] = rowval; colptr[end] = N+1
    nzval  = Vector{Float64}(N); fill!(nzval, 1/N)
    return Operator(SparseMatrixCSC(N,N,colptr,rowval,nzval),(N,))
end
