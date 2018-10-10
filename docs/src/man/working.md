```@meta
DocTestSetup = quote
    using Schrodinger, Compat.SparseArrays
end
```

# Working with States and Operators

Now that we know how to create states and operators, we would like to know what we can do with them. This section describes the fundamental quantum mechanical operations that can be performed in Schrodinger.jl.

## Acting on States with Operators

 - Multiplying operators and states
 - Applying the raising and lowering operator, especially show the lowering operator on the zero state

Operators cannot act directly on density matrices, because operators and density matrices both are linear operators. Instead, we can create a superoperator. These are just operators that act on operators. Mathematically, since density matrices are themselves elements of a vector space (not the same one as for kets, however), we can perform linear operations on them. This is what superoperators do.

## Norms
- Calculating the norm/trace
- Normalization and renormalization

## Tensor Products

## The Partial Trace

## Expectation Values

One of the most fundamental operation in quantum mechanics is the measurement of a physical Hermitian operator, like particle number, or of a non-Hermitian operator. Both operations, though they have different physical implications (physical observables must be Hermitian) are done in the same way: by calculating the [expectation value](https://en.wikipedia.org/wiki/Expectation_value_(quantum_mechanics)) of the operator with respect to a particular state.

Calculating expectation values is done with the [`expect`](@ref) function, both with kets and density matrices, although the familiar mathematical notation can be used as well.

Let's start with a qubit as a simple example.
