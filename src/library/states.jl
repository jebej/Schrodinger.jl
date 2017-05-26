"""
    basis(N, n, dense=false)

Generate a basis state (a.k.a. Fock or number state) ket \$|n⟩\$, in a Hilbert space of size `N`. Note that the size of the Hilbert space must be *at least* `n+1`. The function `fock` is an alias for `basis`.

# Example
```jldoctest
julia> ψ = basis(3,2)
3-d Schrodinger.Ket{SparseVector{Float64,Int64},1} with space dimensions 3:
1.00∠0°|2⟩
```
"""
function basis(N::Integer, n::Integer, dense::Bool=false)
    N>n || throw(ArgumentError("basis level $n is too large for an $N-d space"))
    # Julia is 1-indexed!
    return ket(SparseVector(N,[n+1],[1.0]),(N,),dense)
end

const fock = basis

"""
    coherent(N, α, dense=false; analytic=false)

Generate a coherent state ket \$|α⟩\$, in a Hilbert space of size `N`. To create a coherent density matrix, use the `density` function: `density(coherent(N,n))`.

Two methods can be used for generating a coherent state: via application of a displacment operator on a ground state (the default), or analytically, with the formula

```math
|α⟩ = e^{-\\frac{|α|^2}{2}} \\sum_{n=0}^{N-1} \\frac{α^n}{\\sqrt{n!}} |n⟩.
```

While the operator method will return a normalized ket, the analytic method will not. Both methods converge as `N` gets larger. The analytic method is also much faster, especially for large `N`.

# Example
```jldoctest
julia> coherent(6,0.4+1im)
6-d Schrodinger.Ket{Array{Complex{Float64},1},1} with space dimensions 6:
0.60∠68°|1⟩ + 0.56∠0°|0⟩ + 0.46∠136°|2⟩ + 0.29∠-155°|3⟩ + 0.15∠-87°|4⟩
```
"""
function coherent(N::Integer, α::Number, dense::Bool=true; analytic::Bool=false)
    if analytic
        a = exp(-abs2(α)/2)
        x = [a*α^n/sqrfact(n) for n = 0:N-1]
        return ket(x,(N,),dense)
    else
        return displacementop(N,α,dense)*basis(N,0,dense)
    end
end

"""
    thermal(N, n, dense=false)

Generate a thermal state density matrix \$ρ_n\$ with particle number `n`, in a Hilbert space of size `N`. A thermal state \$ρ_n\$ is a probabilistic mixture of basis states such that the expectation value of the number operator \$\\hat{n}\$ is `n`. Note that this is true only if \$N≫n\$. The returned density matrix is always normalized.

# Example
```jldoctest
julia> N=5;n=0.2;

julia> ρ = thermal(N,n)
5×5 Schrodinger.Density{SparseMatrixCSC{Float64,Int64},1} with space dimensions 5:
 0.833441  0.0       0.0        0.0         0.0
 0.0       0.138907  0.0        0.0         0.0
 0.0       0.0       0.0231511  0.0         0.0
 0.0       0.0       0.0        0.00385852  0.0
 0.0       0.0       0.0        0.0         0.000643087

julia> expect(numberop(N),ρ)
0.19935691318327978
```
"""
function thermal(N::Integer, n::Real, dense::Bool=false)
    β = log(1.0/n+1.0)
    rowval = collect(1:N)
    colptr = Vector{Int}(N+1); colptr[1:N] = rowval; colptr[end] = N+1
    nzval  = normalize!([exp(-β*k) for k = 0:N-1],1)
    return density(SparseMatrixCSC(N,N,colptr,rowval,nzval), (N,), dense)
end

"""
    maxmixed(N, dense=false)

Generate a maximally mixed density matrix. The maximally mixed state is a mixture of basis states with uniform probability.

# Example
```jldoctest
julia> maxmixed(4)
4×4 Schrodinger.Density{SparseMatrixCSC{Float64,Int64},1} with space dimensions 4:
 0.25  0.0   0.0   0.0
 0.0   0.25  0.0   0.0
 0.0   0.0   0.25  0.0
 0.0   0.0   0.0   0.25
```
"""
function maxmixed(N::Integer, dense::Bool=false)
    rowval = collect(1:N)
    colptr = Vector{Int}(N+1); colptr[1:N] = rowval; colptr[end] = N+1
    nzval  = Vector{Float64}(N); fill!(nzval, 1/N)
    return density(SparseMatrixCSC(N,N,colptr,rowval,nzval), (N,), dense)
end
