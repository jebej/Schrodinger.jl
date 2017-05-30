"""
    qzero(N, dims=(N,))

Generate a zero operator for a Hilbert space of size `N`. It is possible to specify the subspace dimensions with the `dims` argument. Returns a sparse matrix.

# Example
```jldoctest
julia> qzero(4,(2,2))
4×4 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},2} with space dimensions 2⊗2:
 0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0
```
"""
function qzero(N::Integer, dims::SDims=(N,))
    rowval = Vector{Int}(0)
    colptr = ones(Int,N+1)
    nzval  = Vector{Float64}(0)
    return Operator(SparseMatrixCSC(N,N,colptr,rowval,nzval), dims)
end


"""
    qeye(N, dims=(N,))

Generate an identity operator for a Hilbert space of size `N`. It is possible to specify the subspace dimensions with the `dims` argument. Returns a sparse matrix.

# Example
```jldoctest
julia> qeye(4,(2,2))
4×4 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},2} with space dimensions 2⊗2:
 1.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0
 0.0  0.0  1.0  0.0
 0.0  0.0  0.0  1.0
```
"""
function qeye(N::Integer, dims::SDims=(N,))
    rowval = collect(1:N)
    colptr = Vector{Int}(N+1); colptr[1:N] = rowval; colptr[end] = N+1
    nzval  = ones(N)
    return Operator(SparseMatrixCSC(N,N,colptr,rowval,nzval), dims)
end

"""
    destroy(N)

Generate a quantum harmonic oscillator lowering (annihilation) operator \$\\hat{a}\$ in a truncated Hilbert space of size `N`. Returns a sparse matrix.

# Example
```jldoctest
julia> destroy(4)
4×4 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with space dimensions 4:
 0.0  1.0  0.0      0.0
 0.0  0.0  1.41421  0.0
 0.0  0.0  0.0      1.73205
 0.0  0.0  0.0      0.0
```
"""
function destroy(N::Integer)
    I = 1:N-1; J = 2:N
    V = [sqrt(i) for i in I]
    return Operator(sparse(I,J,V,N,N), (N,))
end

"""
    create(N)

Generate a quantum harmonic oscillator raising (creation) operator \$\\hat{a}^†\$ in a truncated Hilbert space of size `N`. Returns a sparse matrix.

# Example
```jldoctest
julia> create(4)
4×4 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with space dimensions 4:
 0.0  0.0      0.0      0.0
 1.0  0.0      0.0      0.0
 0.0  1.41421  0.0      0.0
 0.0  0.0      1.73205  0.0
```
"""
function create(N::Integer)
    I = 2:N; J = 1:N-1
    V = [sqrt(j) for j in J]
    return Operator(sparse(I,J,V,N,N), (N,))
end

"""
    numberop(N)

Generate a number operator \$\\hat{n}\$ in a Hilbert space of size `N`. Returns a sparse matrix.

# Example
```jldoctest
julia> numberop(4)
4×4 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with space dimensions 4:
 0.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0
 0.0  0.0  2.0  0.0
 0.0  0.0  0.0  3.0
```
"""
function numberop(N::Integer)
    # "nzval" includes a structural 0 for the [1,1] entry
    rowval = collect(1:N)
    colptr = Vector{Int}(N+1); colptr[1:N] = rowval; colptr[end] = N+1
    nzval  = [float(n) for n = 0:N-1]
    return Operator(SparseMatrixCSC(N,N,colptr,rowval,nzval), (N,))
end

"""
    displacementop(N, α)

Generate a quantum harmonic oscillator displacement operator \$\\hat{D}(α)\$ in a truncated Hilbert space of size `N`. Returns a dense matrix.

```math
\\hat{D}(α) = \\exp\\left(α\\hat{a}^† - α^*\\hat{a}\\right)
```

# Example
```jldoctest
julia> displacementop(3,0.5im)
3×3 Schrodinger.Operator{Array{Complex{Float64},2},1} with space dimensions 3:
   0.88262+0.0im            0.0+0.439802im  -0.166001+0.0im
       0.0+0.439802im  0.647859+0.0im             0.0+0.621974im
 -0.166001+0.0im            0.0+0.621974im    0.76524+0.0im
```
"""
function displacementop(N::Integer, α::Number)
    a = full(data(destroy(N)))
    return Operator(expm(α*a' - α'*a), (N,))
end

"""
    squeezeop(N, z)

Generate a quantum harmonic oscillator squeeze operator \$\\hat{S}(z)\$ in a truncated Hilbert space of size `N`. Returns a dense matrix.

```math
\\hat{S}(z) = \\exp\\left(\\frac{1}{2}\\left(z^*\\hat{a}^2 - z\\hat{a}^{†2}\\right)\\right)
```

# Example
```jldoctest
julia> squeezeop(3,0.5im)
3×3 Schrodinger.Operator{Array{Complex{Float64},2},1} with space dimensions 3:
 0.938148+0.0im       0.0+0.0im       0.0-0.346234im
      0.0+0.0im       1.0+0.0im       0.0+0.0im
      0.0-0.346234im  0.0+0.0im  0.938148+0.0im
```
"""
function squeezeop(N::Integer, z::Number)
    a = full(data(destroy(N)))
    return Operator(expm(0.5*(z'*a^2 - z*a'^2)), (N,))
end
