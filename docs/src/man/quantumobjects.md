```@meta
DocTestSetup  = quote
    using Schrodinger
end
```

# Quantum Objects

This section is an introduction to the basic objects used in Schrodinger.jl. It contains an overview of the quantum object generation capabilities offered, and basic mathematical operations. There are 4 types of basic quantum objects: [`Ket`](@ref), [`Bra`](@ref), [`Density`](@ref), and [`Operator`](@ref). All of these objects can be created from a Julia vector or matrix, as appropriate, with a generating function, or by composition of previously made objects.

Since these objects are very similar to regular vectors and matrices, simple mathematics and functions behave as would be expected.

Functions for creating states and operators are listed in the API sections [State Library](@ref) and [Operator Library](@ref)

## [Creating States] (@id creating_states)

There are three ways of representing quantum states: with `Ket` or `Bra` vectors for pure states, and with `Density` matrices for both pure and mixed states.

We already saw that it is possible to create a pure ket state from a Julia vector using the [`ket`](@ref) function. Kets (and bras) are by default stored as sparse vectors. Schrodinger.jl exposes a few functions to generate common states. These functions are listed in the table below; click on the function name for more details.

| Function           | Type      | Notes                                                               |
|--------------------|-----------|---------------------------------------------------------------------|
| [`basis`](@ref)    | `Ket`     | A simple basis vector. The function `fock` is an alias to this one. |
| [`coherent`](@ref) | `Ket`     | A quantum harmonic oscillator coherent state.                       |
| [`maxmixed`](@ref) | `Density` | The maximally mixed state.                                          |
| [`thermal`](@ref)  | `Density` | A thermal state.                                                    |

### Kets

A simple computational state like $$|3⟩$$ can be created with the [`basis`](@ref) function. [`basis`](@ref) takes two arguments: the dimension of the Hilbert space, and the level. Remember that the ground state is given by the zeroth level.

Let us create a `Ket` for a three-level atom in the first excited `e1` state, which is level "1".

```jldoctest threelevel
julia> e1 = basis(3,1)
3-d Schrodinger.Ket{SparseVector{Float64,Int64},1} with space dimensions 3:
1.00∠0°|1⟩
```

A quantum harmonic oscillator can be in what is called a [coherent state](https://en.wikipedia.org/wiki/Coherent_states). Schrodinger.jl provides a function to create such a state. A coherent state is parameterized by $$α$$, which is a complex number determining the amplitude and phase of the state. Remember that a quantum harmonic oscillator is infinite-dimensional. The state space must therefore be truncated to a finite number of level. The [`coherent`](@ref) function takes two arguments, the truncated space size `N`, and `α`.

```jldoctest coherexample
julia> α = 1.5+1im;

julia> Φ = coherent(10,α)
10-d Schrodinger.Ket{Array{Complex{Float64},1},1} with space dimensions 10:
0.47∠101°|3⟩ + 0.45∠67°|2⟩ + 0.42∠135°|4⟩ + 0.35∠34°|1⟩ + 0.34∠168°|5⟩ +...
```

A coherent state is a superposition of number states, which is evident when displayed in the number basis. Note the three dots at the end of the line: Schrodinger.jl only displays the 5 largest components of a `Ket` vector. You can see the full vector with the [`full`](@ref) function:

```jldoctest coherexample
julia> full(Φ)
10-element Array{Complex{Float64},1}:
 0.196912+1.38778e-16im
          0.295368+0.196912im
          0.174045+0.417707im
        -0.0904438+0.462268im
         -0.298846+0.301357im
         -0.335859+0.0686346im
         -0.231869-0.094343im
        -0.0987909-0.145533im
      -0.000832521-0.0994845im
         0.0472052-0.0721079im
```

A *mixed* state is a probabilistic mixture of *pure* states, and it is important to understand the difference between the two. For example, we can create a superposition between two state of a three-level atom by adding kets together:

```jldoctest threelevel
julia> ψ = e1 + basis(3,0)
3-d Schrodinger.Ket{SparseVector{Float64,Int64},1} with space dimensions 3:
1.00∠0°|0⟩ + 1.00∠0°|1⟩
```

!!! note
    Notice that the coefficients of the new state `ψ` add up to 2. By default, Schrodinger.jl does not renormalize states. This is because adding, for example, three states together would incur two renormalization steps (one per addition) and the resulting state would most likely not be what was desired. Instead, you must add up the states you want in the desired proportions, and then use the [`normalize!`](@ref) function.

Let's make sure that this state is normalized:

```jldoctest threelevel
julia> normalize!(ψ)
3-d Schrodinger.Ket{SparseVector{Float64,Int64},1} with space dimensions 3:
0.71∠0°|0⟩ + 0.71∠0°|1⟩
```

This pure state now represents a physical quantum superposition.

### Density Matrices

Let's now imagine that we have a device that creates three-level atoms, but every time you ask for an one, the machine creates an atom in the state `ψ` with probability 1/3, and in the state `e1` with probability 2/3. After pressing the "new atom" button and obtaining a fresh atom, but *before* looking at it, the state of that atom is unknown. This situation describes a mixed state, and such a state can only be represented by a density matrix.

In bra-ket notation, the a pure state can be transformed in a density matrix by multiplying it on the right with its dual bra: $$|ψ⟩⟨ψ|$$. This is done with the complex transpose operation in Schrodinger.jl. This is therefore how we create the correct state for our mystery atom:

```jldoctest threelevel
julia> ρ = 1/3 * ψ*ψ' + 2/3 * e1*e1'
3×3 Schrodinger.Density{SparseMatrixCSC{Float64,Int64},1} with space dimensions 3:
 0.166667  0.166667  0.0
 0.166667  0.833333  0.0
 0.0       0.0       0.0
```

Notice that because the probabilities 1/2 and 2/3 add up to 1, the density matrix is already properly normalized: its trace is one. If that had not been the case, we could have normalized the density matrix with the `normalize!` function again.

Density matrices can also be created directly from a matrix or from a ket with the [`density`](@ref) function:

```jldoctest threelevel
julia> ρ += density(basis(3,2))
3×3 Schrodinger.Density{SparseMatrixCSC{Float64,Int64},1} with space dimensions 3:
 0.166667  0.166667  0.0
 0.166667  0.833333  0.0
 0.0       0.0       1.0

julia> normalize!(ρ)
3×3 Schrodinger.Density{SparseMatrixCSC{Float64,Int64},1} with space dimensions 3:
 0.0833333  0.0833333  0.0
 0.0833333  0.416667   0.0
 0.0        0.0        0.5
```


## Creating Operators

Operators are used to act on quantum states, either continuously, through time evolution under a Hamiltonian, or discretely. As mentioned previously, kets are element of a Hilbert space. Operators are not elements of that space, they *act* on elements to take them to other elements.

Operators can be created from a Julia matrix with the [`operator`](@ref) function and are by default stored as sparse matrices. As with states, Schrodinger.jl contains functions to create common operators:

| Function                 | Type       | Notes                                                   |
|--------------------------|------------|---------------------------------------------------------|
| [`qzero`](@ref)          | `Operator` | The zero operator.                                      |
| [`qeye`](@ref)           | `Operator` | The identity operator.                                  |
| [`numberop`](@ref)       | `Operator` | The particle number operator.                           |
| [`destroy`](@ref)        | `Operator` | The quantum harmonic oscillator lowering operator.      |
| [`create`](@ref)         | `Operator` | The quantum harmonic oscillator raising operator.       |
| [`displacementop`](@ref) | `Operator` | The quantum harmonic oscillator displacement operator.  |
| [`squeezeop`](@ref)      | `Operator` | The quantum harmonic oscillator squeeze operator.       |

Schordinger.jl also exposes the 3 Pauli matrices, the identity operator, and the raising and lowering operators for two-level systems (qubits) as built-in constants. Those are `σx`, `σy`, `σz`, `σ0`, `σ₊`, and `σ₋`. Note that unlike QuTiP, the qubit raising operator will raise $$|0⟩$$ to $$|1⟩$$.

New operators can be constructed from existing ones by adding them or multiplying them together or with numbers.

```jldoctest
julia> a = destroy(5)
5×5 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with space dimensions 5:
 0.0  1.0  0.0      0.0      0.0
 0.0  0.0  1.41421  0.0      0.0
 0.0  0.0  0.0      1.73205  0.0
 0.0  0.0  0.0      0.0      2.0
 0.0  0.0  0.0      0.0      0.0

julia> a'*a + 1/2
5×5 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with space dimensions 5:
 0.5  0.0  0.0  0.0  0.0
 0.0  1.5  0.0  0.0  0.0
 0.0  0.0  2.5  0.0  0.0
 0.0  0.0  0.0  3.5  0.0
 0.0  0.0  0.0  0.0  4.5

julia> a' + a
5×5 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with space dimensions 5:
 0.0  1.0      0.0      0.0      0.0
 1.0  0.0      1.41421  0.0      0.0
 0.0  1.41421  0.0      1.73205  0.0
 0.0  0.0      1.73205  0.0      2.0
 0.0  0.0      0.0      2.0      0.0
```

!!! note
    Adding and substracting numbers to and from operators adds (substracts) the identity matrix multiplied by that number.

## Basic Mathematical Operations

Basic mathematics with kets, density matrices and operators in Schrodinger.jl is very similar to regular linear algebra with vectors and matrices. This is to be expected, kets are elements of a Hilbert space, which is a vector space, and operators are transformations that take kets to kets. The only difference is that when performing operations between quantum objects, their subspace dimensions must be identical. This condition is explained in more details in the next section.

### Algebra

All basic algebra functions work as expected:

```jldoctest
julia> basis(2,0) + basis(2,1)
2-d Schrodinger.Ket{SparseVector{Float64,Int64},1} with space dimensions 2:
1.00∠0°|0⟩ + 1.00∠0°|1⟩

julia> basis(3,0) + 1
3-d Schrodinger.Ket{SparseVector{Float64,Int64},1} with space dimensions 3:
2.00∠0°|0⟩ + 1.00∠0°|1⟩ + 1.00∠0°|2⟩

julia> 2.5im*basis(2,0)
2-d Schrodinger.Ket{SparseVector{Complex{Float64},Int64},1} with space dimensions 2:
2.50∠90°|0⟩

julia> thermal(4,0.3)/2 + density(coherent(4,1))/2
4×4 Schrodinger.Density{Array{Float64,2},1} with space dimensions 4:
 0.569363   0.184874   0.124977   0.0911074
 0.184874   0.275112   0.125807   0.0917127
 0.124977   0.125807   0.105588   0.0619989
 0.0911074  0.0917127  0.0619989  0.049937

julia> numberop(4) + 1/2
4×4 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with space dimensions 4:
 0.5  0.0  0.0  0.0
 0.0  1.5  0.0  0.0
 0.0  0.0  2.5  0.0
 0.0  0.0  0.0  3.5

julia> create(4)^2
4×4 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with space dimensions 4:
 0.0      0.0      0.0  0.0
 0.0      0.0      0.0  0.0
 1.41421  0.0      0.0  0.0
 0.0      2.44949  0.0  0.0
```

!!! note
    As explained previously, adding and substracting numbers to and from density matrices and operators adds (substracts) the identity matrix multiplied by that number. Adding and substracting quantum objects might also lead to non-normalized states. See the [Norms](@ref) section for more details.

### Functions

Many other mathematical functions are available and work as expected:

 - `exp`, `sqrt`, `log`
 - trig functions are missing for now!
 - `real`, `imag`, `abs`, `abs2`
 - `ctranspose`, `conj`, `transpose`

## Other Quantum Objects

There exist other quantum objects, like [`Liouvillian`](@ref)s and [`Propagator`](@ref)s, but those will be discussed in later sections.
