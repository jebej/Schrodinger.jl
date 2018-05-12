var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Schrodinger.jl-Documentation-1",
    "page": "Home",
    "title": "Schrodinger.jl Documentation",
    "category": "section",
    "text": "Schrodinger.jl is a package for quantum simulations. Currently it is focused on time dynamics, but will hopefully be expanded to perform all kind of operations."
},

{
    "location": "index.html#Introduction-1",
    "page": "Home",
    "title": "Introduction",
    "category": "section",
    "text": "This website serves as documentation for the Schrodinger.jl package. If you are just getting started with Schrodinger.jl, it is recommended that you first read the manual (or at least the Getting Started page). The different sections of the manual are easily accessible on the left sidebar.A few tutorials will soon be available.If you are looking for a particular function, you may browse the API by topic or use the search functionality."
},

{
    "location": "index.html#Contributing-1",
    "page": "Home",
    "title": "Contributing",
    "category": "section",
    "text": "Contributions are highly welcomed! Please see the repository on Github and feel free to open issues if you find a bug or if you feel a feature is missing."
},

{
    "location": "index.html#Author-1",
    "page": "Home",
    "title": "Author",
    "category": "section",
    "text": "This package was written by Jérémy Béjanin. If you find it useful, drop me a line!"
},

{
    "location": "man/gettingstarted.html#",
    "page": "Getting Started",
    "title": "Getting Started",
    "category": "page",
    "text": "DocTestSetup  = quote\n    using Schrodinger\nend"
},

{
    "location": "man/gettingstarted.html#Getting-Started-1",
    "page": "Getting Started",
    "title": "Getting Started",
    "category": "section",
    "text": ""
},

{
    "location": "man/gettingstarted.html#First-Steps-1",
    "page": "Getting Started",
    "title": "First Steps",
    "category": "section",
    "text": "The first step to using Schrodinger.jl is to install it. This is easy to do with Julia\'s package manager. Type Pkg.clone(\"https://github.com/jebej/Schrodinger.jl.git\") at the Julia REPL to download the package to your computer. From there, the package can be used or imported like any other Julia package:using Schrodinger"
},

{
    "location": "man/gettingstarted.html#Quantum-States-1",
    "page": "Getting Started",
    "title": "Quantum States",
    "category": "section",
    "text": "In quantum mechanics, one of the most fundamental object for describing the state of a system is a state vector. State vectors represent pure quantum states, as opposed to mixed quantum stated, but we will get to that later.Schrodinger.jl uses the ubiquitous \"bra-ket\" formalism to describe pure states. A ket is nothing more than a normal (column) vector, which in linear algebra lingo are elements of a vector space. In quantum mechanics, the vector space within which ket vectors live is called a Hilbert space.Kets (and their dual, bras) are therefore finite-, or infinite-dimensional vectors. To create a ket in Schrodinger.jl, use the Ket function with a vector as an argument:julia> g = Ket([1,0])\n2-d Schrodinger.Ket{Array{Float64,1},1} with dimensions 2\n1.00∠0°|0⟩The output is a bit busy, so let us go through it.The input to the Ket function was a 2-d vector [1,0]. This can be seen of the first line of the output, which starts with \"2-d\". The next term is the type of the object, which is a Schrodinger (the package name) Ket (the type itself). The Ket type is parameterized by two values, which are seen within the curly brackets and separated by a comma: first, the type of the underlying data, which here is a 1-d Array of Float64 values (by default, kets and other objects are stored in the format you give to the constructor, although the elements will be converted to floating point values), and second, the total number of subspaces in the full Hilbert space that the Ket lives in. Here, the number is simply 1. The end of the line states the subspace dimensions, but since we have only 1 subspace with dimension 2, this is not very interesting.The second line prints the vector in bra-ket polar notation. Since the vector we passed had the entry \"1\" in the \"zeroth\" (which is the first) dimension, the ket we get is simply 0.Schrodinger.jl only supports finite dimensional Hilbert spaces. If the physical system you want to describe is infinite-dimensional, it will need to be truncated.note: Note\nAll objects in Schrodinger.jl are expressed in the computational or number basis. This means that the ground state is 0, and exited states are numbered starting from 1: 1 2 3To learn more about quantum states in Schrodinger.jl, including mixed states represented by density matrices, please see the next section."
},

{
    "location": "man/gettingstarted.html#Operators-1",
    "page": "Getting Started",
    "title": "Operators",
    "category": "section",
    "text": "States are useful, but we need to do something with them. This is what an Operator is for. Note that Operators have the same Julia type as density matrices (the Operator type), but they can be non-Hermitian, and are in general not normalized.Operators act on elements of a Hilbert space (that is, on kets) to modify them. An operator is thus a like a function that takes as input a ket, and returns a new one. The natural representation for an operator is a matrix, but in Schrodinger.jl you need to use the Operator type, which stores a matrix and other important information about the operator.Arbitrary operators can of course be created, but let\'s take a look at one that is built-in, the _x operator:julia> σx\n2×2 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with dimensions 2\n 0.0  1.0\n 1.0  0.0Notice that the first line of the output is very similar to that of the ket we created above. It lists the dimensions of the matrix, the type (which lists the type of the underlying matrix and the number of subspaces), and the space dimensions (which again is just a single 2-d space).The state g that we created in the previous section is a ground state with the same dimensions. Thus, the _x operator can act on it! This is done simply by multiplying the two objects, with the operator acting to the right on the ket:julia> σx*g\n2-d Schrodinger.Ket{Array{Float64,1},1} with dimensions 2\n1.00∠0°|1⟩As expected, the output is a Ket, but notice the state is now 1! By acting on the ground state 0 with the _x operator, we obtained the excited state. This is because the _x operator is the \"flip\" operator. It takes 0 to 1, and 1 to 0. If we apply _x twice then, we get 0 back:julia> σx*σx*g\n2-d Schrodinger.Ket{Array{Float64,1},1} with dimensions 2\n1.00∠0°|0⟩"
},

{
    "location": "man/gettingstarted.html#Simple-Dynamics-1",
    "page": "Getting Started",
    "title": "Simple Dynamics",
    "category": "section",
    "text": "Now that we have a state and an operator, we can perform some time dynamics! The _x operator, in the context of a spin-1/2 system, can represent a transverse magnetic field. In such a situation, a particle starting in the ground state will undergo sinusoidal oscillations between 0 and 1 due to the action of the field. Let\'s simulate it!We first set up the Hamiltonian, assuming our field has an angular frequency =102 (i.e. 1 Hz). If we look at a 2 sec timespan, we should thus see 2 full periods. To measure the value of the spin at each instant in time, we choose the -_z operator as our observable. The minus sign ensures that the 0 state is the lowest energy one (again, because we are in the computational basis).ω = 1.0*2π # angular frequency\nH = ω/2*σx # Hamiltonian\nt = (0.0,2.0) # timespan\nO = -σz # observable\n# output\n2×2 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with dimensions 2\n -1.0  0.0\n  0.0  1.0note: Note\nSchrodinger.jl uses units where  is equal to 1. Make sure that your Hamiltonian is expressed in units of angular frequency, not energy. If you do have a Hamiltonian expressed in energy units, just use the scale! function: scale!(H,1/ħ). The variable ħ is exported by the module and so can be used as-is.We can now pass all three arguments (H, g and O) to the sesolve function (Schrodinger Equation solver) to solve for the time dynamics! We also pass a keyword argument saveat to make sure we have enough points. As can be seen, the results match with theory:res = sesolve(H, g, t, [O], saveat=linspace(0,2,101))\nreal.(res.evals) ≈ -cos.(ω.*res.times) # check against theory\n# output\ntrueLet\'s plot the results!!isdir(\"img\") && mkdir(\"img\")\nusing Schrodinger, PyPlot\nres = sesolve(π*σx, basis(2,0), (0.0,2.0), [-σz], saveat=2/200)\nfigure(figsize=(8,4.5), dpi=100);\nplot(res.times,real.(res.evals)); xlabel(\"time (s)\"); legend([\"\\$⟨-σ_z⟩\\$\"]); grid()\ntight_layout(true); savefig(joinpath(\"img\",\"gettingstarted-plot.svg\"))using PyPlot\nplot(res.times,real.(res.evals)); xlabel(\"time (s)\"); legend([\"\\$⟨-σ_z⟩\\$\"]); grid()(Image: spin-1/2 system oscillations)As we predicted, the system oscillates between -1, the expectation value of -_z when in the ground state, and 1 when in the excited state.This concludes the first section of the manual. Hopefully you now know enough to get started with simple quantum operations. If you would like to learn more about the other features of Schrodinger.jl, keep reading!"
},

{
    "location": "man/quantumobjects.html#",
    "page": "Quantum Objects",
    "title": "Quantum Objects",
    "category": "page",
    "text": "DocTestSetup  = quote\n    using Schrodinger\nend"
},

{
    "location": "man/quantumobjects.html#Quantum-Objects-1",
    "page": "Quantum Objects",
    "title": "Quantum Objects",
    "category": "section",
    "text": "This section is an introduction to the basic objects used in Schrodinger.jl. It contains an overview of the quantum object generation capabilities offered, and basic mathematical operations. There are 3 types of basic quantum objects: Ket, Bra, and Operator. All of these objects can be created from a Julia vector or matrix, as appropriate, with a generating function, or by composition of previously made objects.Since these objects are very similar to regular vectors and matrices, simple mathematics and functions behave as would be expected.Functions for creating states and operators are listed in the API sections State Library and Operator Library"
},

{
    "location": "man/quantumobjects.html#creating_states-1",
    "page": "Quantum Objects",
    "title": "Creating States",
    "category": "section",
    "text": "There are three ways of representing quantum states: with Ket or Bra vectors for pure states, and with density Operators for both pure and mixed states. Note that a density operator is just a normalized, Hermitian Operator.We already saw that it is possible to create a pure ket state from a Julia vector using the Ket function. Kets (and bras) are by default stored as sparse vectors. Schrodinger.jl exposes a few functions to generate common states. These functions are listed in the table below; click on the function name for more details.Function Type Notes\nbasis sparse Ket A simple basis vector. The function fock is an alias to this one.\ncoherent dense Ket A quantum harmonic oscillator coherent state.\nmaxmixed sparse Operator The maximally mixed state.\nthermal sparse Operator A thermal state."
},

{
    "location": "man/quantumobjects.html#Kets-1",
    "page": "Quantum Objects",
    "title": "Kets",
    "category": "section",
    "text": "A simple computational state like 3 can be created with the basis function. basis takes two arguments: the dimension of the Hilbert space, and the level. Remember that the ground state is given by the zeroth level.Let us create a Ket for a three-level atom in the first excited e1 state, which is level \"1\".julia> e1 = basis(3,1)\n3-d Schrodinger.Ket{SparseVector{Float64,Int64},1} with dimensions 3\n1.00∠0°|1⟩A quantum harmonic oscillator can be in what is called a coherent state. Schrodinger.jl provides a function to create such a state. A coherent state is parameterized by , which is a complex number determining the amplitude and phase of the state. Remember that a quantum harmonic oscillator is infinite-dimensional. The state space must therefore be truncated to a finite number of level. The coherent function takes two arguments, the truncated space size N, and α.julia> α = 1.5+1im;\n\njulia> Φ = coherent(10,α)\n10-d Schrodinger.Ket{Array{Complex{Float64},1},1} with dimensions 10\n0.47∠101°|3⟩ + 0.45∠67°|2⟩ + 0.42∠135°|4⟩ + 0.35∠34°|1⟩ + 0.34∠168°|5⟩ +...A coherent state is a superposition of number states, which is evident when displayed in the number basis. Note the three dots at the end of the line: Schrodinger.jl only displays the 5 largest components of a Ket vector. You can see the full vector with the full function:julia> full(Φ)\n10-element Array{Complex{Float64},1}:\n     0.196912+1.38778e-16im\n     0.295368+0.196912im\n     0.174045+0.417707im\n   -0.0904438+0.462268im\n    -0.298846+0.301357im\n    -0.335859+0.0686346im\n    -0.231869-0.094343im\n   -0.0987909-0.145533im\n -0.000832521-0.0994845im\n    0.0472052-0.0721079imA mixed state is a probabilistic mixture of pure states, and it is important to understand the difference between the two. For example, we can create a superposition between two state of a three-level atom by adding kets together:julia> ψ = e1 + basis(3,0)\n3-d Schrodinger.Ket{SparseVector{Float64,Int64},1} with dimensions 3\n1.00∠0°|0⟩ + 1.00∠0°|1⟩note: Note\nNotice that the coefficients of the new state ψ add up to 2. By default, Schrodinger.jl does not renormalize states. This is because adding, for example, three states together would incur two renormalization steps (one per addition) and the resulting state would most likely not be what was desired. Instead, you must add up the states you want in the desired proportions, and then use the normalize! function.Let\'s make sure that this state is normalized:julia> normalize!(ψ)\n3-d Schrodinger.Ket{SparseVector{Float64,Int64},1} with dimensions 3\n0.71∠0°|0⟩ + 0.71∠0°|1⟩This pure state now represents a physical quantum superposition."
},

{
    "location": "man/quantumobjects.html#Density-Matrices-1",
    "page": "Quantum Objects",
    "title": "Density Matrices",
    "category": "section",
    "text": "Let\'s now imagine that we have a device that creates three-level atoms, but every time you ask for an one, the machine creates an atom in the state ψ with probability 1/3, and in the state e1 with probability 2/3. After pressing the \"new atom\" button and obtaining a fresh atom, but before looking at it, the state of that atom is unknown. This situation describes a mixed state, and such a state can only be represented by a density matrix (or density operator).In bra-ket notation, the a pure state can be transformed in a density matrix by multiplying it on the right with its dual bra: . This is done with the complex transpose operation in Schrodinger.jl. This is therefore how we create the correct state for our mystery atom:julia> ρ = 1/3 * ψ*ψ\' + 2/3 * e1*e1\'\n3×3 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with dimensions 3\n 0.166667  0.166667  0.0\n 0.166667  0.833333  0.0\n 0.0       0.0       0.0Notice that because the probabilities 1/2 and 2/3 add up to 1, the matrix is already properly normalized: its trace is one. If that had not been the case, we could have normalized the density operator with the normalize! function again. The density operator of the atom is Hermitian, as it should be. Schrodinger.jl uses the same type (Operator) to represent both density matrices and \"regular\" linear operators.Density matrices can be created directly from a matrix or from a ket with the Operator function:julia> ρ += Operator(basis(3,2))\n3×3 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with dimensions 3\n 0.166667  0.166667  0.0\n 0.166667  0.833333  0.0\n 0.0       0.0       1.0\n\njulia> normalize!(ρ)\n3×3 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with dimensions 3\n 0.0833333  0.0833333  0.0\n 0.0833333  0.416667   0.0\n 0.0        0.0        0.5"
},

{
    "location": "man/quantumobjects.html#Creating-Operators-1",
    "page": "Quantum Objects",
    "title": "Creating Operators",
    "category": "section",
    "text": "Linear operators are used to act on quantum states, either continuously, through time evolution under a Hamiltonian, or discretely. As mentioned previously, kets are element of a Hilbert space. Operators are not elements of that space, they act on elements to take them to other elements.Operators can be created from a Julia matrix with the Operator function and are by default stored as sparse matrices. As with states, Schrodinger.jl contains functions to create common operators:Function Type Notes\nqzero sparse Operator The zero operator.\nqeye sparse Operator The identity operator.\nnumberop sparse Operator The particle number operator.\ndestroy sparse Operator The quantum harmonic oscillator lowering operator.\ncreate sparse Operator The quantum harmonic oscillator raising operator.\ndisplacementop dense Operator The quantum harmonic oscillator displacement operator.\nsqueezeop dense Operator The quantum harmonic oscillator squeeze operator.Schordinger.jl also exposes the 3 Pauli matrices, the identity operator, and the raising and lowering operators for two-level systems (qubits) as built-in constants. Those are σx, σy, σz, σ0, σ₊, and σ₋. Note that unlike QuTiP, the qubit raising operator will raise 0 to 1.New operators can be constructed from existing ones by adding them or multiplying them together or with numbers. Operators can be non-Hermitian, unlike density matrices.julia> a = destroy(5)\n5×5 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with dimensions 5\n 0.0  1.0  0.0      0.0      0.0\n 0.0  0.0  1.41421  0.0      0.0\n 0.0  0.0  0.0      1.73205  0.0\n 0.0  0.0  0.0      0.0      2.0\n 0.0  0.0  0.0      0.0      0.0\n\njulia> a\'*a + 1/2\n5×5 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with dimensions 5\n 0.5  0.0  0.0  0.0  0.0\n 0.0  1.5  0.0  0.0  0.0\n 0.0  0.0  2.5  0.0  0.0\n 0.0  0.0  0.0  3.5  0.0\n 0.0  0.0  0.0  0.0  4.5\n\njulia> a\' + a\n5×5 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with dimensions 5\n 0.0  1.0      0.0      0.0      0.0\n 1.0  0.0      1.41421  0.0      0.0\n 0.0  1.41421  0.0      1.73205  0.0\n 0.0  0.0      1.73205  0.0      2.0\n 0.0  0.0      0.0      2.0      0.0note: Note\nAdding and substracting numbers to and from operators adds (substracts) the identity matrix multiplied by that number."
},

{
    "location": "man/quantumobjects.html#Basic-Mathematical-Operations-1",
    "page": "Quantum Objects",
    "title": "Basic Mathematical Operations",
    "category": "section",
    "text": "Basic mathematics with kets, density matrices and operators in Schrodinger.jl is very similar to regular linear algebra with vectors and matrices. This is to be expected, kets are elements of a Hilbert space, which is a vector space, and operators are transformations that take kets to kets. The only difference is that when performing operations between quantum objects, their subspace dimensions must be identical. This condition is explained in more details in the next section."
},

{
    "location": "man/quantumobjects.html#Algebra-1",
    "page": "Quantum Objects",
    "title": "Algebra",
    "category": "section",
    "text": "All basic algebra functions work as expected:julia> basis(2,0) + basis(2,1)\n2-d Schrodinger.Ket{SparseVector{Float64,Int64},1} with dimensions 2\n1.00∠0°|0⟩ + 1.00∠0°|1⟩\n\njulia> basis(3,0) + 1\n3-d Schrodinger.Ket{SparseVector{Float64,Int64},1} with dimensions 3\n2.00∠0°|0⟩ + 1.00∠0°|1⟩ + 1.00∠0°|2⟩\n\njulia> 2.5im*basis(2,0)\n2-d Schrodinger.Ket{SparseVector{Complex{Float64},Int64},1} with dimensions 2\n2.50∠90°|0⟩\n\njulia> thermal(4,0.3)/2 + Operator(coherent(4,1))/2\n4×4 Schrodinger.Operator{Array{Float64,2},1} with dimensions 4\n 0.569363   0.184874   0.124977   0.0911074\n 0.184874   0.275112   0.125807   0.0917127\n 0.124977   0.125807   0.105588   0.0619989\n 0.0911074  0.0917127  0.0619989  0.049937\n\njulia> numberop(4) + 1/2\n4×4 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with dimensions 4\n 0.5  0.0  0.0  0.0\n 0.0  1.5  0.0  0.0\n 0.0  0.0  2.5  0.0\n 0.0  0.0  0.0  3.5\n\njulia> create(4)^2\n4×4 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with dimensions 4\n 0.0      0.0      0.0  0.0\n 0.0      0.0      0.0  0.0\n 1.41421  0.0      0.0  0.0\n 0.0      2.44949  0.0  0.0note: Note\nAs explained previously, adding and substracting numbers to and from operators adds (substracts) the identity matrix multiplied by that number. Adding and substracting quantum objects might also lead to non-normalized states. See the Norms section for more details."
},

{
    "location": "man/quantumobjects.html#Functions-1",
    "page": "Quantum Objects",
    "title": "Functions",
    "category": "section",
    "text": "Many other mathematical functions are available and work as expected:exp, sqrt, log\ntrig functions are missing for now!\nreal, imag, abs, abs2\nctranspose, conj, transpose"
},

{
    "location": "man/quantumobjects.html#Other-Quantum-Objects-1",
    "page": "Quantum Objects",
    "title": "Other Quantum Objects",
    "category": "section",
    "text": "There exist other quantum objects, like Liouvillians and Propagators, but those will be discussed in later sections."
},

{
    "location": "man/working.html#",
    "page": "Working with States and Operators",
    "title": "Working with States and Operators",
    "category": "page",
    "text": "DocTestSetup  = quote\n    using Schrodinger\nend"
},

{
    "location": "man/working.html#Working-with-States-and-Operators-1",
    "page": "Working with States and Operators",
    "title": "Working with States and Operators",
    "category": "section",
    "text": "Now that we know how to create states and operators, we would like to know what we can do with them. This section describes the fundamental quantum mechanical operations that can be performed in Schrodinger.jl."
},

{
    "location": "man/working.html#Acting-on-States-with-Operators-1",
    "page": "Working with States and Operators",
    "title": "Acting on States with Operators",
    "category": "section",
    "text": "Multiplying operators and states\nApplying the raising and lowering operator, especially show the lowering operator on the zero stateOperators cannot act directly on density matrices, because operators and density matrices both are linear operators. Instead, we can create a superoperator. These are just operators that act on operators. Mathematically, since density matrices are themselves elements of a vector space (not the same one as for kets, however), we can perform linear operations on them. This is what superoperators do."
},

{
    "location": "man/working.html#Norms-1",
    "page": "Working with States and Operators",
    "title": "Norms",
    "category": "section",
    "text": "Calculating the norm/trace\nNormalization and renormalization"
},

{
    "location": "man/working.html#Tensor-Products-1",
    "page": "Working with States and Operators",
    "title": "Tensor Products",
    "category": "section",
    "text": ""
},

{
    "location": "man/working.html#The-Partial-Trace-1",
    "page": "Working with States and Operators",
    "title": "The Partial Trace",
    "category": "section",
    "text": ""
},

{
    "location": "man/working.html#Expectation-Values-1",
    "page": "Working with States and Operators",
    "title": "Expectation Values",
    "category": "section",
    "text": "One of the most fundamental operation in quantum mechanics is the measurement of a physical Hermitian operator, like particle number, or of a non-Hermitian operator. Both operations, though they have different physical implications (physical observables must be Hermitian) are done in the same way: by calculating the expectation value of the operator with respect to a particular state.Calculating expectation values is done with the expect function, both with kets and density matrices, although the familiar mathematical notation can be used as well.Let\'s start with a qubit as a simple example."
},

{
    "location": "man/dynamics.html#",
    "page": "Time Evolution and Dynamics",
    "title": "Time Evolution and Dynamics",
    "category": "page",
    "text": ""
},

{
    "location": "man/dynamics.html#Time-Evolution-and-Dynamics-1",
    "page": "Time Evolution and Dynamics",
    "title": "Time Evolution and Dynamics",
    "category": "section",
    "text": ""
},

{
    "location": "examples/DRAG.html#",
    "page": "DRAG",
    "title": "DRAG",
    "category": "page",
    "text": "DocTestSetup  = quote\n    using Schrodinger\nend!isdir(\"img\") && mkdir(\"img\")"
},

{
    "location": "examples/DRAG.html#DRAG-1",
    "page": "DRAG",
    "title": "DRAG",
    "category": "section",
    "text": "This example shows how to implement a NOT gate on a qubit using a Gaussian pulse. We then extend this method to a three-level slightly anharmonic energy spectrum with nearest level coupling. As we will see, the gate error increases due to leakage into the third level. To remedy this, we implement Derivative Removal by Adiabatic Gate (DRAG) which offers better gate fidelity than an ordinary Gaussian pulse (see the References)."
},

{
    "location": "examples/DRAG.html#A-simple-NOT-Gate-1",
    "page": "DRAG",
    "title": "A simple NOT Gate",
    "category": "section",
    "text": "First, we will apply a simple NOT gate to a qubit in the ground state. The Hamiltonian for our qubit in the lab frame can be written as (11+(t)_x) where  is the transition energy, _x is the Pauli-X operator, and (t)=^x(t)cos(_dt) represents our control of the system using a drive frequency _d. Any control ^x(t) such that the integral of ^x over the total gate time equals  will result in a complete inversion. It is common to use a Gaussian shaped π-pulse to implement a NOT gate. For our purposes, we will find it more convenient to work in the rotating frame with respect to our drive frequency _d. When this frequency is resonant with the qubit frequency , the Hamiltonian is given by ^x(t)_x2. If we are going to solve the time dynamics for this system, we also have to define our pulse in the rotating frame.using Schrodinger\n\nfunction rotgaussianpulse(t::Real,p::Vector)\n    # normalized pulse centered on t=0, begins and ends at 0\n    tg  = p[1] # gate time\n    σ   = p[2] # standard deviation\n    A   = p[3] # amplitude\n    B   = inv(√(2π)*σ*erf(tg/(√(8)*σ))-tg*gaussian(0.5tg,σ)) # normalize\n    Ɛˣ = A*B*(gaussian(t,σ)-gaussian(0.5tg,σ))\n    return Ɛˣ\nend\n\nH = σx/2 # same as create(2)/2 + destroy(2)/2; using natural units where ħ=1We are now ready to solve for the time dynamics and plot our results. We\'ll use parameters that are typical for Transmon qubits.tg = 6e-9 # gate time 6ns\nσ = 0.5tg # standard deviation of gaussian pulse\ng = basis(2,0) # begin in ground state\ntspan = (-tg/2,tg/2) # pulse centered at t=0\n\nres1 = sesolve([qzero(2),(H,rotgaussianpulse,[tg,σ,π])],g,tspan,saveat=tg/200)\n\nusing PyPlot\nfigure(figsize=(8,4.5), dpi=100);\nplot(res1.times*1e9,levelprobs(res1.states)); grid();\nxlabel(\"Time (ns)\"); ylabel(\"Level Probabilities\");\nlegend([\"Ground State\", \"Excited State\"]);\ntight_layout(true); savefig(joinpath(\"img\",\"qubitNOT.svg\")); # hide(Image: qubit NOT gate)As is expected, the system moves from the ground state to the excited state."
},

{
    "location": "examples/DRAG.html#Three-level-system-1",
    "page": "DRAG",
    "title": "Three-level system",
    "category": "section",
    "text": "In reality, one must worry about leakage into other states of the system, especially when dealing with short gate times. Let\'s see how our Gaussian pulse performs on a three-level anharmonic system. We will again work in the rotating frame at resonance with the qubit frequency. This system will have an anharmonicity , which is the detuning of the 2nd excited state with respect to the drive frequency, and an additional parameter  describing the relative strength of the 1-2 transition compared to the 0-1 transition (see [1] for more details). Let\'s create the Hamiltonian and perform the same time evolution.Δ = 2π*(-400e6) # anharmonicity\nλ = √2 # relative transition strength\nΠ₂ = basis(3,2) * basis(3,2)\' # projector for the 2nd level\nHc = Δ*Π₂ # constant Hamiltonian\nHd = create(3)/2 + destroy(3)/2 # affected by Ɛˣ\ng = basis(3,0)\n\nres2 = sesolve([Hc,(Hd,rotgaussianpulse,[tg,σ,π])],g,tspan,saveat=tg/200)\nfigure(figsize=(8,4.5), dpi=100);\nplot(res2.times*1e9,levelprobs(res2.states)); grid()\nxlabel(\"Time (ns)\"); ylabel(\"Level Probabilities\");\nlegend([\"Ground State\", \"1st Excited State\", \"2nd Excited State\"]);\ntight_layout(true); savefig(joinpath(\"img\",\"3levelNOT.svg\")); # hide(Image: 3-level NOT gate)Instead of working perfectly, our system leaks into the 2nd energy level (we\'ll quantify this later). This becomes problematic when we\'re trying perform useful computations. DRAG is one possible remedy where introducing drive detuning and a second quadrature control, that is (t)=^x(t)cos(_dt)+^y(t)sin(_dt), helps eliminate some of the leakage. Let\'s see how much of an improvement we get. We\'ll need to define a few more functions for the new controls.function DRAGx(t::Real,p::Vector)\n    # Ɛˣ for fifth order DRAG\n    tg  = p[1] # gate time\n    σ   = p[2] # standard deviation\n    A   = p[3] # amplitude\n    Δ   = p[4] # anharmonicity\n    λ   = p[5] # relative strength of transitions\n    Ɛˣ = rotgaussianpulse(t,[tg,σ,A]) + (λ^2-4)*rotgaussianpulse(t,[tg,σ,A])^3/(8*Δ^2) - (13λ^4-76λ^2+112)*rotgaussianpulse(t,[tg,σ,A])^5/(128Δ^4)\n    return Ɛˣ\nend\n\nfunction DRAGy(t::Real,p::Vector)\n    # Ɛʸ for fifth order DRAG\n    tg  = p[1] # gate time\n    σ   = p[2] # standard deviation\n    A   = p[3] # amplitude\n    Δ   = p[4] # anharmonicity\n    λ   = p[5] # relative strength of transitions\n    B   = inv(√(2π)*σ*erf(tg/(√(8)*σ))-tg*gaussian(0.5tg,σ))\n    Ɛˣ′ = A*B*Schrodinger.gaussianprime(t,σ)\n    Ɛʸ = -Ɛˣ′/Δ + 33*(λ^2-2)*rotgaussianpulse(t,[tg,σ,A])^2*Ɛˣ′/(24*Δ^3)\n    return Ɛʸ\nend\n\nfunction dydet(t::Real,p::Vector)\n    # dynamical detuning for fifth order DRAG\n    tg  = p[1] # gate time\n    σ   = p[2] # standard deviation\n    A   = p[3] # amplitude\n    Δ   = p[4] # anharmonicity\n    λ   = p[5] # relative strength of transitions\n    δ₁ = (λ^2-4)*rotgaussianpulse(t,[tg,σ,A])^2/(4*Δ) - (λ^4-7λ^2+12)*rotgaussianpulse(t,[tg,σ,A])^4/(16Δ^3)\n    return δ₁\nend\n\nΠ₁ = basis(3,1) * basis(3,1)\' # projector for 1st level\nHdet = Π₁ # affected by dynamical detuning\nHdx = create(3)/2 + destroy(3)/2 # affected by Ɛˣ\nHdy = im*create(3)/2 - im*destroy(3)/2 # affected by Ɛʸ\n\nres3 = sesolve([Hc,(Hdet,dydet,[tg,σ,π,Δ,λ]),(Hdx,DRAGx,[tg,σ,π,Δ,λ]),(Hdy,DRAGy,[tg,σ,π,Δ,λ])],g,tspan,saveat=tg/200)\nfigure(figsize=(8,4.5), dpi=100);\nplot(res3.times*1e9,levelprobs(res3.states)); grid()\nxlabel(\"Time (ns)\"); ylabel(\"Level Probabilities\");\nlegend([\"Ground State\", \"1st Excited State\", \"2nd Excited State\"]);\ntight_layout(true); savefig(joinpath(\"img\",\"3levelDRAG.svg\")); # hide(Image: 3-level DRAG)There is a noticeable improvement over the simple Gaussian pulse, but how much better is the new gate?"
},

{
    "location": "examples/DRAG.html#Fidelity-1",
    "page": "DRAG",
    "title": "Fidelity",
    "category": "section",
    "text": "We can calculate the fidelity of our gates by comparing their output to the ideal case. Our gates behave ideally when =0 and there is no leakage into the 2nd excited state. Using the same measure of error as in [1], we can take the overall gate fidelity to be the average of gate fidelities when using the 6 axial states on the Bloch sphere as inputs.tgs = (2:0.5:10)*1E-9 # gate times\nFg_res = Matrix{Float64}(length(tgs),2) # initialize matrix for number of solves\naxialkets  = [normalize!(Ket([1,1,0])),   # +X\n              normalize!(Ket([1,-1,0])),  # -X\n              normalize!(Ket([1,im,0])),  # +Y\n              normalize!(Ket([1,-im,0])), # -Y\n              normalize!(Ket([1,0,0])),   # +Z\n              normalize!(Ket([0,1,0]))]   # -Z\naxialops = Operator.(axialkets) # need density operators too\n\nUideal = qeye(3); Uideal[1:2,1:2] = data(σx)\nfor (i,tg) in enumerate(tgs)\n    sum1 = 0 # sum of Gaussian gate fidelities\n    sum2 = 0 # sum of DRAG gate fidelities\n    for j = 1:6\n        ket = axialkets[j]\n        op = axialops[j]\n        tspan = (-tg/2,tg/2)\n        # Gaussian\n        res4_1 = sesolve([Hc,(Hd,rotgaussianpulse,[tg,σ,π])],ket,tspan)\n        sum1 += trace(Uideal*op*Uideal\'*(res4_1.states[end]*res4_1.states[end]\'))\n        # DRAG\n        res4_2 = sesolve([Hc,(Hdet,dydet,[tg,σ,π,Δ,λ]),(Hdx,DRAGx,[tg,σ,π,Δ,λ]),(Hdy,DRAGy,[tg,σ,π,Δ,λ])],ket,tspan)\n        sum2 += trace(Uideal*op*Uideal\'*(res4_2.states[end]*res4_2.states[end]\'))\n    end\n    Fg_res[i,:] = [sum1/6 sum2/6] # take average\nend\n\nfigure(figsize=(8,4.5), dpi=100);\nplot(tgs*1E9,1.-Fg_res); ylim([10E-8,1]); grid()\ntitle(\"Average gate fidelity averaging over all input states\");\nyscale(\"log\"); xlabel(\"Gate Time (ns)\"); ylabel(\"Gate Error 1-Fg\");\nlegend([\"Gaussian\",\"DRAG 5th Order\"]);\ntight_layout(true); savefig(joinpath(\"img\",\"fidelities.svg\")); # hide(Image: NOT-gate fidelities)For a gate time of 6ns, taking advantage of DRAG results in a gate error that is 2 orders of magnitude less than when using Gaussian pulses."
},

{
    "location": "examples/DRAG.html#References-1",
    "page": "DRAG",
    "title": "References",
    "category": "section",
    "text": "Motzoi, F., Gambetta, J. M., Rebentrost, P., & Wilhelm, F. K. (2009). Simple Pulses for Elimination of Leakage in Weakly Nonlinear Qubits. Physical Review Letters, 103(11), 110501. 10.1103/PhysRevLett.103.110501 or arXiv:0901.0534"
},

{
    "location": "api/quobj.html#",
    "page": "Quantum Object Types",
    "title": "Quantum Object Types",
    "category": "page",
    "text": "DocTestSetup  = quote\n    using Schrodinger\nend"
},

{
    "location": "api/quobj.html#Schrodinger.Bra",
    "page": "Quantum Object Types",
    "title": "Schrodinger.Bra",
    "category": "type",
    "text": "Bra(x, dims=(length(x),))\n\nBra vector type. The dual vector to the Ket.\n\nThe Bra type has two fields, data and dims, which store the vector data and the subspace dimensions. A Bra, like a Ket or an Operator is parameterized by the number of subspaces it lives in. Two different kets must have the same system dimensions in order to be added together.\n\nIt is possible to normalize the bra vector after construction with the normalize! function.\n\n\n\n"
},

{
    "location": "api/quobj.html#Schrodinger.Ket",
    "page": "Quantum Object Types",
    "title": "Schrodinger.Ket",
    "category": "type",
    "text": "Ket(x, dims=(length(x),))\n\nConstruct a ket state vector from the vector x. A vector of length N will by default be assumed to be an element of a single Hilbert space of dimension N. If the vector is an element of a tensor product of Hilbert spaces, the dimensions can be defined manually by passing a tuple of subspace dimensions dims. In that case, prod(dims) must equal length(x). By default, the vector is stored in sparse format.\n\nThe Ket type has two fields, data and dims, which store the vector data and the subspace dimensions. A Ket, like a Bra or an Operator is parameterized by the number of subspaces it lives in. Two different kets must have the same system dimensions in order to be added together.\n\nIt is possible to normalize the ket vector after construction with the normalize! function.\n\nExample\n\njulia> ψ = normalize!(Ket([1,1]))\n2-d Schrodinger.Ket{Array{Float64,1},1} with dimensions 2\n0.71∠0°|0⟩ + 0.71∠0°|1⟩\n\n\n\n"
},

{
    "location": "api/quobj.html#Schrodinger.Operator",
    "page": "Quantum Object Types",
    "title": "Schrodinger.Operator",
    "category": "type",
    "text": "Operator(B, dims=(size(B,1),))\n\nConstruct a linear operator from the matrix B. An N×N matrix will by default be assumed to describe an operator that acts on a single Hilbert space of dimension N. If the matrix represents a linear operator on a tensor product of Hilbert spaces, the dimensions can be defined manually by passing a tuple of subspace dimensions dims. In that case, prod(dims) must equal size(B,1).\n\nThe Operator type has two fields, data and dims, which store the matrix data and the subspace dimensions. An Operator, like a Ket or a Bra, is parameterized by the number of subspaces it lives in. Two different density matrices must have the same system dimensions in order to be added together. An Operator may or may not be Hermitian.\n\nExample\n\njulia> σ = Operator([0 -im ; im 0])\n2×2 Schrodinger.Operator{Array{Complex{Float64},2},1} with dimensions 2\n 0.0+0.0im  0.0-1.0im\n 0.0+1.0im  0.0+0.0im\n\n\n\n"
},

{
    "location": "api/quobj.html#Quantum-Object-Types-1",
    "page": "Quantum Object Types",
    "title": "Quantum Object Types",
    "category": "section",
    "text": "Modules = [Schrodinger]\nPages   = [\"types.jl\"]\nOrder   = [:type, :function]\nPrivate = false"
},

{
    "location": "api/states.html#",
    "page": "State Library",
    "title": "State Library",
    "category": "page",
    "text": "DocTestSetup  = quote\n    using Schrodinger\nend"
},

{
    "location": "api/states.html#Schrodinger.basis-Tuple{Integer,Tuple{Vararg{Int64,D}} where D}",
    "page": "State Library",
    "title": "Schrodinger.basis",
    "category": "method",
    "text": "basis(N, n)\n\nGenerate a basis state (a.k.a. Fock or number state) ket n, in a Hilbert space of size N. Note that the size of the Hilbert space must be at least n+1. The function fock is an alias for basis.\n\nReturns a sparse vector.\n\nExample\n\njulia> ψ = basis(3,2)\n3-d Schrodinger.Ket{SparseVector{Float64,Int64},1} with dimensions 3\n1.00∠0°|2⟩\n\n\n\n"
},

{
    "location": "api/states.html#Schrodinger.coherent",
    "page": "State Library",
    "title": "Schrodinger.coherent",
    "category": "function",
    "text": "coherent(N, α, analytic=false)\n\nGenerate a coherent state ket , in a Hilbert space of size N. To create a coherent density operator, use the Operator function: Operator(coherent(N,n)).\n\nTwo methods can be used for generating a coherent state: via application of a displacment operator on a ground state (the default), or analytically, with the formula\n\n = e^-frac^22 sum_n=0^N-1 frac^nsqrtn n\n\nWhile the operator method will return a normalized ket, the analytic method will not. Both methods converge as N gets larger. The analytic method is also much faster, especially for large N.\n\nReturns a dense vector.\n\nExample\n\njulia> coherent(6,0.4+1im)\n6-d Schrodinger.Ket{Array{Complex{Float64},1},1} with dimensions 6\n0.60∠68°|1⟩ + 0.56∠0°|0⟩ + 0.46∠136°|2⟩ + 0.29∠-155°|3⟩ + 0.15∠-87°|4⟩\n\n\n\n"
},

{
    "location": "api/states.html#Schrodinger.ket-Union{Tuple{N}, Tuple{Tuple{Vararg{Int64,N}},Int64}, Tuple{Tuple{Vararg{Int64,N}}}} where N",
    "page": "State Library",
    "title": "Schrodinger.ket",
    "category": "method",
    "text": "ket(state,dims=2)\n\nGenerate a state ket from a tuple of basis levels and a tuple of corresponding space dimensions. Note that each space dimension must be larger than the level by at least 1. If only an integer is passed to dims, all basis levels will have the same dimension.\n\nReturns a sparse vector.\n\nExample\n\njulia> ψ = ket((3,0,1),(5,2,3))\n30-d Schrodinger.Ket{SparseVector{Float64,Int64},3} with dimensions 5⊗2⊗3\n1.00∠0°|3,0,1⟩\n\nSee also: qb, for qubit states.\n\n\n\n"
},

{
    "location": "api/states.html#Schrodinger.maxentangled",
    "page": "State Library",
    "title": "Schrodinger.maxentangled",
    "category": "function",
    "text": "maxentangled(n,N=2)\n\nGenerate a maximally entangled state between n N-d systems:\n\nphi=sum_j=0^N-1frac1sqrtNj^n\n\nTracing out all but one of the entangled systems results in a maximally mixed state.\n\nExample\n\njulia> ψ = maxentangled(3,4)\n64-d Schrodinger.Ket{SparseVector{Float64,Int64},3} with dimensions 4⊗4⊗4\n0.50∠0°|0,0,0⟩ + 0.50∠0°|1,1,1⟩ + 0.50∠0°|2,2,2⟩ + 0.50∠0°|3,3,3⟩\n\njulia> ptrace(ψ,(1,3))\n4×4 Schrodinger.Operator{Array{Float64,2},1} with dimensions 4\n 0.25  0.0   0.0   0.0\n 0.0   0.25  0.0   0.0\n 0.0   0.0   0.25  0.0\n 0.0   0.0   0.0   0.25\n\n\n\n"
},

{
    "location": "api/states.html#Schrodinger.maxmixed-Tuple{Integer}",
    "page": "State Library",
    "title": "Schrodinger.maxmixed",
    "category": "method",
    "text": "maxmixed(N)\n\nGenerate a maximally mixed density matrix. The maximally mixed state is a mixture of basis states with uniform probability.\n\nReturns a sparse matrix.\n\nExample\n\njulia> maxmixed(4)\n4×4 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with dimensions 4\n 0.25  0.0   0.0   0.0\n 0.0   0.25  0.0   0.0\n 0.0   0.0   0.25  0.0\n 0.0   0.0   0.0   0.25\n\n\n\n"
},

{
    "location": "api/states.html#Schrodinger.qb-Tuple{Vararg{Int64,N} where N}",
    "page": "State Library",
    "title": "Schrodinger.qb",
    "category": "method",
    "text": "qb(q1,q2,q3...)\n\nGenerate a qubit state from the given argument list. This function is similar to ket, except that the state levels are passed with separate arguments instead of a tuple.\n\nReturns a sparse vector.\n\nExample\n\njulia> Ψ⁻ = normalize!(qb(0,1) - qb(1,0))\n4-d Schrodinger.Ket{SparseVector{Float64,Int64},2} with dimensions 2⊗2\n0.71∠0°|0,1⟩ + 0.71∠180°|1,0⟩\n\n\n\n"
},

{
    "location": "api/states.html#Schrodinger.thermal-Tuple{Integer,Real}",
    "page": "State Library",
    "title": "Schrodinger.thermal",
    "category": "method",
    "text": "thermal(N, n)\n\nGenerate a thermal state density matrix _n with particle number n, in a Hilbert space of size N. A thermal state _n is a probabilistic mixture of basis states such that the expectation value of the number operator hatn is n. Note that this is true only if Nn. The returned density matrix is always normalized.\n\nReturns a sparse matrix.\n\nExample\n\njulia> N=5; n=0.2;\n\njulia> ρ = thermal(N,n)\n5×5 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with dimensions 5\n 0.833441  0.0       0.0        0.0         0.0\n 0.0       0.138907  0.0        0.0         0.0\n 0.0       0.0       0.0231511  0.0         0.0\n 0.0       0.0       0.0        0.00385852  0.0\n 0.0       0.0       0.0        0.0         0.000643087\n\njulia> expect(numberop(N),ρ)\n0.19935691318327978\n\n\n\n"
},

{
    "location": "api/states.html#State-Library-1",
    "page": "State Library",
    "title": "State Library",
    "category": "section",
    "text": "Modules = [Schrodinger]\nPages   = [\"states.jl\"]\nOrder   = [:function]\nPrivate = false"
},

{
    "location": "api/operators.html#",
    "page": "Operator Library",
    "title": "Operator Library",
    "category": "page",
    "text": "DocTestSetup = quote\n    using Schrodinger\nend"
},

{
    "location": "api/operators.html#Schrodinger.create-Tuple{Integer}",
    "page": "Operator Library",
    "title": "Schrodinger.create",
    "category": "method",
    "text": "create(N)\n\nGenerate a quantum harmonic oscillator raising (creation) operator hata^ in a truncated Hilbert space of size N. Returns a sparse matrix.\n\nExample\n\njulia> create(4)\n4×4 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with dimensions 4\n 0.0  0.0      0.0      0.0\n 1.0  0.0      0.0      0.0\n 0.0  1.41421  0.0      0.0\n 0.0  0.0      1.73205  0.0\n\n\n\n"
},

{
    "location": "api/operators.html#Schrodinger.destroy-Tuple{Integer}",
    "page": "Operator Library",
    "title": "Schrodinger.destroy",
    "category": "method",
    "text": "destroy(N)\n\nGenerate a quantum harmonic oscillator lowering (annihilation) operator hata in a truncated Hilbert space of size N. Returns a sparse matrix.\n\nExample\n\njulia> destroy(4)\n4×4 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with dimensions 4\n 0.0  1.0  0.0      0.0\n 0.0  0.0  1.41421  0.0\n 0.0  0.0  0.0      1.73205\n 0.0  0.0  0.0      0.0\n\n\n\n"
},

{
    "location": "api/operators.html#Schrodinger.displacementop-Tuple{Integer,Number}",
    "page": "Operator Library",
    "title": "Schrodinger.displacementop",
    "category": "method",
    "text": "displacementop(N, α)\n\nGenerate a quantum harmonic oscillator displacement operator hatD() in a truncated Hilbert space of size N. Returns a dense matrix.\n\nhatD() = expleft(hata^ - ^*hataright)\n\nExample\n\njulia> displacementop(3,0.5im)\n3×3 Schrodinger.Operator{Array{Complex{Float64},2},1} with dimensions 3\n   0.88262+0.0im            0.0+0.439802im  -0.166001+0.0im\n       0.0+0.439802im  0.647859+0.0im             0.0+0.621974im\n -0.166001+0.0im            0.0+0.621974im    0.76524+0.0im\n\n\n\n"
},

{
    "location": "api/operators.html#Schrodinger.numberop-Tuple{Integer}",
    "page": "Operator Library",
    "title": "Schrodinger.numberop",
    "category": "method",
    "text": "numberop(N)\n\nGenerate a number operator hatn in a Hilbert space of size N. Returns a sparse matrix.\n\nExample\n\njulia> numberop(4)\n4×4 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with dimensions 4\n 0.0  0.0  0.0  0.0\n 0.0  1.0  0.0  0.0\n 0.0  0.0  2.0  0.0\n 0.0  0.0  0.0  3.0\n\n\n\n"
},

{
    "location": "api/operators.html#Schrodinger.projectorop-Union{Tuple{Integer,AbstractArray{T,1}}, Tuple{T}} where T<:Integer",
    "page": "Operator Library",
    "title": "Schrodinger.projectorop",
    "category": "method",
    "text": "projectorop(N,S)\n\nGenerate a projector on the subspaces defined by an integer or a vector/range of integers S:\n\nP = sum_iS ii\n\nExample\n\njulia> projectorop(5,[1,3])\n5×5 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with dimensions 5\n 0.0  0.0  0.0  0.0  0.0\n 0.0  1.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  1.0  0.0\n 0.0  0.0  0.0  0.0  0.0\n\n\n\n"
},

{
    "location": "api/operators.html#Schrodinger.qeye",
    "page": "Operator Library",
    "title": "Schrodinger.qeye",
    "category": "function",
    "text": "qeye(N, dims=(N,))\n\nGenerate an identity operator for a Hilbert space of size N. It is possible to specify the subspace dimensions with the dims argument. Returns a sparse matrix.\n\nExample\n\njulia> qeye(4,(2,2))\n4×4 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},2} with dimensions 2⊗2\n 1.0  0.0  0.0  0.0\n 0.0  1.0  0.0  0.0\n 0.0  0.0  1.0  0.0\n 0.0  0.0  0.0  1.0\n\n\n\n"
},

{
    "location": "api/operators.html#Schrodinger.qzero",
    "page": "Operator Library",
    "title": "Schrodinger.qzero",
    "category": "function",
    "text": "qzero(N, dims=(N,))\n\nGenerate a zero operator for a Hilbert space of size N. It is possible to specify the subspace dimensions with the dims argument. Returns a sparse matrix.\n\nExample\n\njulia> qzero(4,(2,2))\n4×4 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},2} with dimensions 2⊗2\n 0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0\n\n\n\n"
},

{
    "location": "api/operators.html#Schrodinger.squeezeop-Tuple{Integer,Number}",
    "page": "Operator Library",
    "title": "Schrodinger.squeezeop",
    "category": "method",
    "text": "squeezeop(N, z)\n\nGenerate a quantum harmonic oscillator squeeze operator hatS(z) in a truncated Hilbert space of size N. Returns a dense matrix.\n\nhatS(z) = expleft(frac12left(z^*hata^2 - zhata^2right)right)\n\nExample\n\njulia> squeezeop(3,0.5im)\n3×3 Schrodinger.Operator{Array{Complex{Float64},2},1} with dimensions 3\n 0.938148+0.0im       0.0+0.0im       0.0-0.346234im\n      0.0+0.0im       1.0+0.0im       0.0+0.0im\n      0.0-0.346234im  0.0+0.0im  0.938148+0.0im\n\n\n\n"
},

{
    "location": "api/operators.html#Schrodinger.rand_unitary",
    "page": "Operator Library",
    "title": "Schrodinger.rand_unitary",
    "category": "function",
    "text": "rand_unitary(N, dims=(N,))\n\nGenerate a Haar distributed random unitary operator for a Hilbert space of size N. It is possible to specify the subspace dimensions with the dims argument. Returns a dense matrix.\n\nExample\n\njulia> U = rand_unitary(4,(2,2));\n\njulia> U\'*U ≈ qeye(4,(2,2))\ntrue\n\n\n\n"
},

{
    "location": "api/operators.html#Operator-Library-1",
    "page": "Operator Library",
    "title": "Operator Library",
    "category": "section",
    "text": "Modules = [Schrodinger]\nPages   = [\"operators.jl\",\"random.jl\"]\nOrder   = [:function]\nPrivate = false"
},

{
    "location": "api/functions.html#",
    "page": "Function Library",
    "title": "Function Library",
    "category": "page",
    "text": "DocTestSetup  = quote\n    using Schrodinger\nend"
},

{
    "location": "api/functions.html#Schrodinger.expect-Tuple{Schrodinger.Operator,Schrodinger.Ket}",
    "page": "Function Library",
    "title": "Schrodinger.expect",
    "category": "method",
    "text": "expect(σ,ψ), expect(σ,ρ)\n\nCompute the expectation value of an operator  given a state ket  or a density matrix . The expectation value is defined as\n\nbeginalign*\nE() =  \nE() = textrmtr()\nendalign*\n\nA specialized method exists for vector of Ket or Operator inputs\n\n\n\n"
},

{
    "location": "api/functions.html#Schrodinger.fidelity-Tuple{Schrodinger.Operator,Schrodinger.Operator}",
    "page": "Function Library",
    "title": "Schrodinger.fidelity",
    "category": "method",
    "text": "fidelity(ρ,σ), fidelity(ρ,ψ), fidelity(ψ,ϕ)\n\nCompute the fidelity between density matrices  and , a density matrix  and a ket , or two kets  and . The fidelity in those three cases is defined as\n\nbeginalign*\nF() = textrmtrsqrt^12^12 \nF() = sqrt \nF() = leftright\nendalign*\n\nSee also fidelity2, which is the square of the state fidelity.\n\n\n\n"
},

{
    "location": "api/functions.html#Schrodinger.fidelity2-Tuple{Schrodinger.Operator,Schrodinger.Operator}",
    "page": "Function Library",
    "title": "Schrodinger.fidelity2",
    "category": "method",
    "text": "fidelity2(ρ,ψ)\n\nCompute the Uhlmann state fidelity between density matrices  and , a density matrix  and a ket , or two kets  and . The Uhlmann state fidelity is simply defined as the square of the \"regular\" state fidelity.\n\n\n\n"
},

{
    "location": "api/functions.html#Schrodinger.levelprobs-Union{Tuple{N}, Tuple{Schrodinger.Ket{T,N}}, Tuple{T}} where N where T",
    "page": "Function Library",
    "title": "Schrodinger.levelprobs",
    "category": "method",
    "text": "levelprobs(ψ), levelprobs(ψ,s)\n\nCompute the level occupation probabilities. For a Ket, this simply corresponds to the absolute square of the amplitude of each level. For an Operator, the function returns the diagonal.\n\nA system index, or vector of indices, can be passed as a second argument. In that case, the full system will first be partial traced to keep only the desired index. Level occupation probabilities are then calculated from the resulting reduced density matrix. If a vector of indices is passed, occupation probabilities are calculated for a fully reduced density matrix for each index.\n\nA specialized method exists for vector of Ket or Operator inputs.\n\n\n\n"
},

{
    "location": "api/functions.html#Schrodinger.ptrace-Tuple{Schrodinger.Ket,Any}",
    "page": "Function Library",
    "title": "Schrodinger.ptrace",
    "category": "method",
    "text": "ptrace(ψ, out)\n\nCompute the partial trace of a state Ket or Bra ψ by tracing out the subsystems specified by out. Returns a density matrix. Multiple subsystems can be traced out by passing a sorted tuple of subsystem indices.\n\nExample\n\njulia> Φ₊ = normalize!(basis(2,0)⊗basis(2,0) + basis(2,1)⊗basis(2,1)) # Bell pair\n4-d Schrodinger.Ket{SparseVector{Float64,Int64},2} with dimensions 2⊗2\n0.71∠0°|0,0⟩ + 0.71∠0°|1,1⟩\n\njulia> ptrace(Φ₊,1) # trace out qubit 1\n2×2 Schrodinger.Operator{Array{Float64,2},1} with dimensions 2\n 0.5  0.0\n 0.0  0.5\n\n\n\n"
},

{
    "location": "api/functions.html#Schrodinger.ptrace-Tuple{Schrodinger.Operator,Any}",
    "page": "Function Library",
    "title": "Schrodinger.ptrace",
    "category": "method",
    "text": "ptrace(ρ, out)\n\nCompute the partial trace of a linear Operator ρ by tracing out the subsystems specified by out. Multiple subsystems can be traced out by passing a sorted tuple of subsystem indices.\n\nExample\n\nΦ₊ = normalize!(basis(2,0)⊗basis(2,0) + basis(2,1)⊗basis(2,1)) # Bell pair\nΨ₊ = normalize!(basis(2,0)⊗basis(2,1) + basis(2,1)⊗basis(2,0)) # Bell pair\nρ  = 0.25 * Operator(Φ₊) + 0.75 * Operator(Ψ₊) # density matrix\nptrace(ρ,2) # trace out qubit 2\n# output\n2×2 Schrodinger.Operator{Array{Float64,2},1} with dimensions 2\n 0.5  0.0\n 0.0  0.5\n\n\n\n"
},

{
    "location": "api/functions.html#Function-Library-1",
    "page": "Function Library",
    "title": "Function Library",
    "category": "section",
    "text": "Modules = [Schrodinger]\nPages   = [\"special.jl\",\"ptrace.jl\"]\nOrder   = [:function]\nPrivate = false"
},

]}
