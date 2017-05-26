"""
    qzero(N, dims=(N,), dense=false)

Generate a zero operator for a Hilbert space of size `N`. It is possible to specify the subspace dimensions with the `dims` argument.

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
function qzero(N::Integer, dims::SchroDims=(N,), dense::Bool=false)
    return operator(spzeros(N,N),dims,dense)
end

"""
    qeye(N, dims=(N,), dense=false)

Generate an identity operator for a Hilbert space of size `N`. It is possible to specify the subspace dimensions with the `dims` argument.

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
function qeye(N::Integer, dims::SchroDims=(N,), dense::Bool=false)
    rowval = collect(1:N)
    colptr = Vector{Int}(N+1); colptr[1:N] = rowval; colptr[end] = N+1
    nzval  = ones(N)
    return operator(SparseMatrixCSC(N,N,colptr,rowval,nzval), dims, dense)
end

"""
    destroy(N, dense=false)

Generate a quantum harmonic oscillator lowering (annihilation) operator \$\\hat{a}\$ in a truncated Hilbert space of size `N`.

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
function destroy(N::Integer, dense::Bool=false)
    I = 1:N-1; J = 2:N
    V = [sqrt(i) for i in I]
    return operator(sparse(I,J,V,N,N), (N,), dense)
end

"""
    create(N, dense=false)

Generate a quantum harmonic oscillator raising (creation) operator \$\\hat{a}^†\$ in a truncated Hilbert space of size `N`.

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
function create(N::Integer, dense::Bool=false)
    I = 2:N; J = 1:N-1
    V = [sqrt(j) for j in J]
    return operator(sparse(I,J,V,N,N), (N,), dense)
end

"""
    numberop(N, dense=false)

Generate a number operator \$\\hat{n}\$ in a Hilbert space of size `N`.

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
function numberop(N::Integer, dense::Bool=false)
    # "nzval" includes a structural 0 for the [1,1] entry
    rowval = collect(1:N)
    colptr = Vector{Int}(N+1); colptr[1:N] = rowval; colptr[end] = N+1
    nzval  = [float(n) for n = 0:N-1]
    return operator(SparseMatrixCSC(N,N,colptr,rowval,nzval),(N,),dense)
end

"""
    displacementop(N, α, dense=true)

Generate a quantum harmonic oscillator displacement operator \$\\hat{D}(α)\$ in a truncated Hilbert space of size `N`.

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
function displacementop(N::Integer, α::Number, dense::Bool=true)
    a = destroy(N,true).data
    return operator(expm(α*a' - α'*a), (N,), dense)
end

"""
    squeezeop(N, z, dense=true)

Generate a quantum harmonic oscillator squeeze operator \$\\hat{S}(z)\$ in a truncated Hilbert space of size `N`.

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
function squeezeop(N::Integer, z::Number, dense::Bool=true)
    a = destroy(N,true).data
    return operator(expm(0.5*(z'*a^2 - z*a'^2)), (N,), dense)
end
