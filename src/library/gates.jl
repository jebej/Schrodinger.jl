module Gate
using ..Schrodinger
using SparseArrays

"""
    Gate.rotation(θ, n=(1,0,0))

Generate a qubit rotation operator giving a `θ` rad rotation about an axis defined by the vector ``\\vec{n}``. Note that `n` will be normalized, allowing for inputs like `(1,1,0)`. The rotation operator is defined as

```math
\\hat{R}_{\\vec{n}}(θ) = \\exp\\left(−iθ\\vec{n}⋅\\vec{σ}/2\\right) = \\cos\\frac{θ}{2}𝕀 - i\\sin\\frac{θ}{2}(n_x σ_x + n_y σ_y + n_z σ_z).
```
"""
function rotation(θ::Real,n::NTuple{3,Real}=(1,0,0))
    a = 1/√(n[1]^2+n[2]^2+n[3]^2)
    nx,ny,nz = a*n[1],a*n[2],a*n[3]
    c = cos(θ/2); s = sin(θ/2)
    R = Matrix{Complex{typeof(nz*s)}}(undef,2,2)
    R[1,1] = complex(c,-nz*s)
    R[2,1] = complex(ny*s,-nx*s)
    R[1,2] = complex(-ny*s,-nx*s)
    R[2,2] = complex(c,nz*s)
    return Operator(R,(2,))
end

"""
    Gate.H

The Hadamard gate:
```math
H = \\frac{1}{\\sqrt 2}
\\begin{pmatrix}
1 & 1 \\\\
1 & -1
\\end{pmatrix}
```

The Hadamard is a ``π`` rotation about the axis ``\\vec{n} = (1,0,1)``, plus a global ``i`` phase.
"""
const H = Operator([√0.5 √0.5; √0.5 -√0.5],(2,))

"""
    Gate.S

The phase gate:
```math
S =
\\begin{pmatrix}
1 & 0 \\\\
0 & i
\\end{pmatrix}
```

The phase gate is the square of the T gate: ``S = T^2``.
"""
const S = Operator(sparse([1 0; 0 1im]),(2,))

"""
    Gate.T

The T gate:
```math
T =
\\begin{pmatrix}
1 & 0 \\\\
0 & \\exp(iπ/4)
\\end{pmatrix}
```

The phase gate is the square of the T gate: ``S = T^2``.
"""
const T = Operator(sparse([1 0; 0 complex(√0.5,√0.5)]),(2,))

const cNOT = Operator(sparse([1 0 0 0;0 1 0 0;0 0 0 1;0 0 1 0]),(2,2))
const rcNOT = Operator(sparse([0 1 0 0;1 0 0 0;0 0 1 0;0 0 0 1]),(2,2))
const cZ = Operator(sparse([1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 -1]),(2,2))

end
