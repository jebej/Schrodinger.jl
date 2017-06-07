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
    "text": "The first step to using Schrodinger.jl is to install it. This is easy to do with Julia's package manager. Type Pkg.clone(\"https://github.com/jebej/Schrodinger.jl.git\") at the Julia REPL to download the package to your computer. From there, the package can be used or imported like any other Julia package:using Schrodinger"
},

{
    "location": "man/gettingstarted.html#Quantum-States-1",
    "page": "Getting Started",
    "title": "Quantum States",
    "category": "section",
    "text": "In quantum mechanics, one of the most fundamental object for describing the state of a system is a state vector. State vectors represent pure quantum states, as opposed to mixed quantum stated, but we will get to that later.Schrodinger.jl uses the ubiquitous \"bra-ket\" formalism to describe pure states. A ket is nothing more than a normal (column) vector, which in linear algebra lingo are elements of a vector space. In quantum mechanics, the vector space within which ket vectors live is called a Hilbert space.Kets (and their dual, bras) are therefore finite-, or infinite-dimensional vectors. To create a ket in Schrodinger.jl, use the Ket function with a vector as an argument:julia> g = Ket([1,0])\n2-d Schrodinger.Ket{Array{Float64,1},1} with space dimensions 2:\n1.00∠0°|0⟩The output is a bit busy, so let us go through it.The input to the Ket function was a 2-d vector [1,0]. This can be seen of the first line of the output, which starts with \"2-d\". The next term is the type of the object, which is a Schrodinger (the package name) Ket (the type itself). The Ket type is parameterized by two values, which are seen within the curly brackets and separated by a comma: first, the type of the underlying data, which here is a 1-d Array of Float64 values (by default, kets and other objects are stored in the format you give to the constructor, although the elements will be converted to floating point values), and second, the total number of subspaces in the full Hilbert space that the Ket lives in. Here, the number is simply 1. The end of the line states the subspace dimensions, but since we have only 1 subspace with dimension 2, this is not very interesting.The second line prints the vector in bra-ket polar notation. Since the vector we passed had the entry \"1\" in the \"zeroth\" (which is the first) dimension, the ket we get is simply 0.Schrodinger.jl only supports finite dimensional Hilbert spaces. If the physical system you want to describe is infinite-dimensional, it will need to be truncated.note: Note\nAll objects in Schrodinger.jl are expressed in the computational or number basis. This means that the ground state is 0, and exited states are numbered starting from 1: 1 2 3To learn more about quantum states in Schrodinger.jl, including mixed states represented by density matrices, please see the next section."
},

{
    "location": "man/gettingstarted.html#Operators-1",
    "page": "Getting Started",
    "title": "Operators",
    "category": "section",
    "text": "States are useful, but we need to do something with them. This is what an Operator is for. Operators act on elements of a Hilbert space (that is, on kets) to modify them. An operator is thus a like a function that takes as input a ket, and returns a new one. The natural representation for an operator is a matrix, but in Schrodinger.jl you need to use the Operator type, which stores a matrix and other important information about the operator.Arbitrary operators can of course be created, but let's take a look at one that is built-in, the _x operator:julia> σx\n2×2 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with space dimensions 2:\n 0.0  1.0\n 1.0  0.0Notice that the first line of the output is very similar to that of the ket we created above. It lists the dimensions of the matrix, the type and the space dimensions (which again is just a single 2-d space).The state g that we created in the previous section is a ground state with the same dimensions. Thus, the _x operator can act on it! This is done simply by multiplying the two objects, with the operator acting to the right on the ket:julia> σx*g\n2-d Schrodinger.Ket{Array{Float64,1},1} with space dimensions 2:\n1.00∠0°|1⟩As expected, the output is a Ket, but notice the state is now 1! By acting on the ground state 0 with the _x operator, we obtained the excited state. This is because the _x operator is the \"flip\" operator. It takes 0 to 1, and 1 to 0. If we apply _x twice then, we get 0 back:julia> σx*σx*g\n2-d Schrodinger.Ket{Array{Float64,1},1} with space dimensions 2:\n1.00∠0°|0⟩"
},

{
    "location": "man/gettingstarted.html#Simple-Dynamics-1",
    "page": "Getting Started",
    "title": "Simple Dynamics",
    "category": "section",
    "text": "Now that we have a state and an operator, we can perform some time dynamics! The _x operator, in the context of a spin-1/2 system, can represent a transverse magnetic field. In such a situation, a particle starting in the ground state will undergo sinusoidal oscillations between 0 and 1 due to the action of the field. Let's simulate it!We first set up the Hamiltonian, assuming our field has an angular frequency =102 (i.e. 1 Hz). If we look at a 2 sec timespan, we should thus see 2 full periods. To measure the value of the spin at each instant in time, we choose the -_z operator as our observable. The minus sign ensures that the 0 state is the lowest energy one (again, because we are in the computational basis).ω = 1.0*2π # angular frequency\nH = ω/2*σx # Hamiltonian\nt = (0.0,2.0) # timespan\nO = -σz # observable\n# output\n2×2 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with space dimensions 2:\n -1.0  0.0\n  0.0  1.0note: Note\nSchrodinger.jl uses units where  is equal to 1. Make sure that your Hamiltonian is expressed in units of angular frequency, not energy. If you do have a Hamiltonian expressed in energy units, just use the scale! function: scale!(H,1/ħ). The variable ħ is exported by the module and so can be used as-is.We can now pass all three arguments (H, g and O) to the sesolve function (Schrodinger Equation solver) to solve for the time dynamics! We also pass a keyword argument saveat to make sure we have enough points. As can be seen, the results match with theory:res = sesolve(H, g, t, [O], saveat=linspace(0,2,101))\nreal.(res.evals) ≈ -cos(ω.*res.times) # check against theory\n# output\ntrueLet's plot the results!print(\"Loading Schrodinger.jl...\")\nusing Schrodinger\nprintln(\" done!\")\nprint(\"Running differential equation solver for the first time...\")\nres = sesolve(π*σx, basis(2,0), (0.0,2.0), [-σz], saveat=linspace(0,2,101))\nprintln(\" done!\")\nprint(\"Loading Plots.jl...\")\nusing Plots\nprintln(\" done!\")\nprint(\"Plotting for the first time...\")\nplot([1,2,3],[1,2,3]);\nprintln(\" done!\")\n!isdir(\"img\") && mkdir(\"img\")using Plots\nplot(res.times, real.(res.evals), xlabel=\"time (s)\", label=\"\\$⟨-σ_z⟩\\$\")\nsavefig(joinpath(\"img\",\"example1-plot.svg\")); nothing # hide(Image: spin-1/2 system oscillations)As we predicted, the system oscillates between -1, the expectation value of -_z when in the ground state, and 1 when in the excited state.This concludes the first section of the manual. Hopefully you now know enough to get started with simple quantum operations. If you would like to learn more about the other features of Schrodinger.jl, keep reading!"
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
    "text": "This section is an introduction to the basic objects used in Schrodinger.jl. It contains an overview of the quantum object generation capabilities offered, and basic mathematical operations. There are 4 types of basic quantum objects: Ket, Bra, Density, and Operator. All of these objects can be created from a Julia vector or matrix, as appropriate, with a generating function, or by composition of previously made objects.Since these objects are very similar to regular vectors and matrices, simple mathematics and functions behave as would be expected.Functions for creating states and operators are listed in the API sections State Library and Operator Library"
},

{
    "location": "man/quantumobjects.html#creating_states-1",
    "page": "Quantum Objects",
    "title": "Creating States",
    "category": "section",
    "text": "There are three ways of representing quantum states: with Ket or Bra vectors for pure states, and with Density matrices for both pure and mixed states.We already saw that it is possible to create a pure ket state from a Julia vector using the Ket function. Kets (and bras) are by default stored as sparse vectors. Schrodinger.jl exposes a few functions to generate common states. These functions are listed in the table below; click on the function name for more details.Function Type Notes\nbasis sparse Ket A simple basis vector. The function fock is an alias to this one.\ncoherent dense Ket A quantum harmonic oscillator coherent state.\nmaxmixed sparse Density The maximally mixed state.\nthermal sparse Density A thermal state."
},

{
    "location": "man/quantumobjects.html#Kets-1",
    "page": "Quantum Objects",
    "title": "Kets",
    "category": "section",
    "text": "A simple computational state like 3 can be created with the basis function. basis takes two arguments: the dimension of the Hilbert space, and the level. Remember that the ground state is given by the zeroth level.Let us create a Ket for a three-level atom in the first excited e1 state, which is level \"1\".julia> e1 = basis(3,1)\n3-d Schrodinger.Ket{SparseVector{Float64,Int64},1} with space dimensions 3:\n1.00∠0°|1⟩A quantum harmonic oscillator can be in what is called a coherent state. Schrodinger.jl provides a function to create such a state. A coherent state is parameterized by , which is a complex number determining the amplitude and phase of the state. Remember that a quantum harmonic oscillator is infinite-dimensional. The state space must therefore be truncated to a finite number of level. The coherent function takes two arguments, the truncated space size N, and α.julia> α = 1.5+1im;\n\njulia> Φ = coherent(10,α)\n10-d Schrodinger.Ket{Array{Complex{Float64},1},1} with space dimensions 10:\n0.47∠101°|3⟩ + 0.45∠67°|2⟩ + 0.42∠135°|4⟩ + 0.35∠34°|1⟩ + 0.34∠168°|5⟩ +...A coherent state is a superposition of number states, which is evident when displayed in the number basis. Note the three dots at the end of the line: Schrodinger.jl only displays the 5 largest components of a Ket vector. You can see the full vector with the full function:julia> full(Φ)\n10-element Array{Complex{Float64},1}:\n 0.196912+1.38778e-16im\n          0.295368+0.196912im\n          0.174045+0.417707im\n        -0.0904438+0.462268im\n         -0.298846+0.301357im\n         -0.335859+0.0686346im\n         -0.231869-0.094343im\n        -0.0987909-0.145533im\n      -0.000832521-0.0994845im\n         0.0472052-0.0721079imA mixed state is a probabilistic mixture of pure states, and it is important to understand the difference between the two. For example, we can create a superposition between two state of a three-level atom by adding kets together:julia> ψ = e1 + basis(3,0)\n3-d Schrodinger.Ket{SparseVector{Float64,Int64},1} with space dimensions 3:\n1.00∠0°|0⟩ + 1.00∠0°|1⟩note: Note\nNotice that the coefficients of the new state ψ add up to 2. By default, Schrodinger.jl does not renormalize states. This is because adding, for example, three states together would incur two renormalization steps (one per addition) and the resulting state would most likely not be what was desired. Instead, you must add up the states you want in the desired proportions, and then use the normalize! function.Let's make sure that this state is normalized:julia> normalize!(ψ)\n3-d Schrodinger.Ket{SparseVector{Float64,Int64},1} with space dimensions 3:\n0.71∠0°|0⟩ + 0.71∠0°|1⟩This pure state now represents a physical quantum superposition."
},

{
    "location": "man/quantumobjects.html#Density-Matrices-1",
    "page": "Quantum Objects",
    "title": "Density Matrices",
    "category": "section",
    "text": "Let's now imagine that we have a device that creates three-level atoms, but every time you ask for an one, the machine creates an atom in the state ψ with probability 1/3, and in the state e1 with probability 2/3. After pressing the \"new atom\" button and obtaining a fresh atom, but before looking at it, the state of that atom is unknown. This situation describes a mixed state, and such a state can only be represented by a density matrix.In bra-ket notation, the a pure state can be transformed in a density matrix by multiplying it on the right with its dual bra: . This is done with the complex transpose operation in Schrodinger.jl. This is therefore how we create the correct state for our mystery atom:julia> ρ = 1/3 * ψ*ψ' + 2/3 * e1*e1'\n3×3 Schrodinger.Density{SparseMatrixCSC{Float64,Int64},1} with space dimensions 3:\n 0.166667  0.166667  0.0\n 0.166667  0.833333  0.0\n 0.0       0.0       0.0Notice that because the probabilities 1/2 and 2/3 add up to 1, the density matrix is already properly normalized: its trace is one. If that had not been the case, we could have normalized the density matrix with the normalize! function again.Density matrices can also be created directly from a matrix or from a ket with the Density function:julia> ρ += Density(basis(3,2))\n3×3 Schrodinger.Density{SparseMatrixCSC{Float64,Int64},1} with space dimensions 3:\n 0.166667  0.166667  0.0\n 0.166667  0.833333  0.0\n 0.0       0.0       1.0\n\njulia> normalize!(ρ)\n3×3 Schrodinger.Density{SparseMatrixCSC{Float64,Int64},1} with space dimensions 3:\n 0.0833333  0.0833333  0.0\n 0.0833333  0.416667   0.0\n 0.0        0.0        0.5"
},

{
    "location": "man/quantumobjects.html#Creating-Operators-1",
    "page": "Quantum Objects",
    "title": "Creating Operators",
    "category": "section",
    "text": "Operators are used to act on quantum states, either continuously, through time evolution under a Hamiltonian, or discretely. As mentioned previously, kets are element of a Hilbert space. Operators are not elements of that space, they act on elements to take them to other elements.Operators can be created from a Julia matrix with the Operator function and are by default stored as sparse matrices. As with states, Schrodinger.jl contains functions to create common operators:Function Type Notes\nqzero sparse Operator The zero operator.\nqeye sparse Operator The identity operator.\nnumberop sparse Operator The particle number operator.\ndestroy sparse Operator The quantum harmonic oscillator lowering operator.\ncreate sparse Operator The quantum harmonic oscillator raising operator.\ndisplacementop dense Operator The quantum harmonic oscillator displacement operator.\nsqueezeop dense Operator The quantum harmonic oscillator squeeze operator.Schordinger.jl also exposes the 3 Pauli matrices, the identity operator, and the raising and lowering operators for two-level systems (qubits) as built-in constants. Those are σx, σy, σz, σ0, σ₊, and σ₋. Note that unlike QuTiP, the qubit raising operator will raise 0 to 1.New operators can be constructed from existing ones by adding them or multiplying them together or with numbers.julia> a = destroy(5)\n5×5 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with space dimensions 5:\n 0.0  1.0  0.0      0.0      0.0\n 0.0  0.0  1.41421  0.0      0.0\n 0.0  0.0  0.0      1.73205  0.0\n 0.0  0.0  0.0      0.0      2.0\n 0.0  0.0  0.0      0.0      0.0\n\njulia> a'*a + 1/2\n5×5 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with space dimensions 5:\n 0.5  0.0  0.0  0.0  0.0\n 0.0  1.5  0.0  0.0  0.0\n 0.0  0.0  2.5  0.0  0.0\n 0.0  0.0  0.0  3.5  0.0\n 0.0  0.0  0.0  0.0  4.5\n\njulia> a' + a\n5×5 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with space dimensions 5:\n 0.0  1.0      0.0      0.0      0.0\n 1.0  0.0      1.41421  0.0      0.0\n 0.0  1.41421  0.0      1.73205  0.0\n 0.0  0.0      1.73205  0.0      2.0\n 0.0  0.0      0.0      2.0      0.0note: Note\nAdding and substracting numbers to and from operators adds (substracts) the identity matrix multiplied by that number."
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
    "text": "All basic algebra functions work as expected:julia> basis(2,0) + basis(2,1)\n2-d Schrodinger.Ket{SparseVector{Float64,Int64},1} with space dimensions 2:\n1.00∠0°|0⟩ + 1.00∠0°|1⟩\n\njulia> basis(3,0) + 1\n3-d Schrodinger.Ket{SparseVector{Float64,Int64},1} with space dimensions 3:\n2.00∠0°|0⟩ + 1.00∠0°|1⟩ + 1.00∠0°|2⟩\n\njulia> 2.5im*basis(2,0)\n2-d Schrodinger.Ket{SparseVector{Complex{Float64},Int64},1} with space dimensions 2:\n2.50∠90°|0⟩\n\njulia> thermal(4,0.3)/2 + Density(coherent(4,1))/2\n4×4 Schrodinger.Density{Array{Float64,2},1} with space dimensions 4:\n 0.569363   0.184874   0.124977   0.0911074\n 0.184874   0.275112   0.125807   0.0917127\n 0.124977   0.125807   0.105588   0.0619989\n 0.0911074  0.0917127  0.0619989  0.049937\n\njulia> numberop(4) + 1/2\n4×4 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with space dimensions 4:\n 0.5  0.0  0.0  0.0\n 0.0  1.5  0.0  0.0\n 0.0  0.0  2.5  0.0\n 0.0  0.0  0.0  3.5\n\njulia> create(4)^2\n4×4 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with space dimensions 4:\n 0.0      0.0      0.0  0.0\n 0.0      0.0      0.0  0.0\n 1.41421  0.0      0.0  0.0\n 0.0      2.44949  0.0  0.0note: Note\nAs explained previously, adding and substracting numbers to and from density matrices and operators adds (substracts) the identity matrix multiplied by that number. Adding and substracting quantum objects might also lead to non-normalized states. See the Norms section for more details."
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
    "text": "One of the most fundamental operation in quantum mechanics is the measurement of a physical Hermitian operator, like particle number, or of a non-Hermitian operator. Both operations, though they have different physical implications (physical observables must be Hermitian) are done in the same way: by calculating the expectation value of the operator with respect to a particular state.Calculating expectation values is done with the expect function, both with kets and density matrices, although the familiar mathematical notation can be used as well.Let's start with a qubit as a simple example."
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
    "category": "Type",
    "text": "Bra(x, dims=(length(x),))\n\nBra vector type. The dual vector to the Ket.\n\nThe Bra type has two fields, data and dims, which store the vector data and the subspace dimensions. A Bra, like a Density matrix or and Operator is parameterized by the number of subspaces it lives in. Two different kets must have the same system dimensions in order to be added together.\n\nIt is possible to normalize the bra vector after construction with the normalize! function.\n\n\n\n"
},

{
    "location": "api/quobj.html#Schrodinger.Density",
    "page": "Quantum Object Types",
    "title": "Schrodinger.Density",
    "category": "Type",
    "text": "Density(A, dims=(size(A,1),))\n\nConstruct a density matrix (a.k.a. density operator) from the Hermitian matrix A. An N×N matrix will by default be assumed to describe a single Hilbert space of dimension N. If the matrix represents a tensor product of Hilbert spaces, the dimensions can be defined manually by passing a tuple of subspace dimensions dims. In that case, prod(dims) must equal size(A,1). By default, the matrix is stored in sparse format.\n\nThe Density type has two fields, data and dims, which store the matrix data and the subspace dimensions. A Density matrix, like a Ket or an Operator, is parameterized by the number of subspaces it lives in. Two different density matrices must have the same system dimensions in order to be added together. A Density matrix is always Hermitian, this is enforced on construction by \"hermitianizing\" the given matrix.\n\nIt is possible to normalize the density matrix after construction with the normalize! function.\n\nExample\n\njulia> A = [1 5 2; 5 1 0; 2 0 2]\n3×3 Array{Int64,2}:\n 1  5  2\n 5  1  0\n 2  0  2\njulia> σ = normalize!(Density(A))\n3×3 Schrodinger.Density{Array{Float64,2},1} with space dimensions 3:\n 0.25  1.25  0.5\n 1.25  0.25  0.0\n 0.5   0.0   0.5\njulia> trace(σ)\n1.0\n\n\n\n"
},

{
    "location": "api/quobj.html#Schrodinger.Ket",
    "page": "Quantum Object Types",
    "title": "Schrodinger.Ket",
    "category": "Type",
    "text": "Ket(x, dims=(length(x),))\n\nConstruct a ket state vector from the vector x. A vector of length N will by default be assumed to be an element of a single Hilbert space of dimension N. If the vector is an element of a tensor product of Hilbert spaces, the dimensions can be defined manually by passing a tuple of subspace dimensions dims. In that case, prod(dims) must equal length(x). By default, the vector is stored in sparse format.\n\nThe Ket type has two fields, data and dims, which store the vector data and the subspace dimensions. A Ket, like a Density matrix or and Operator is parameterized by the number of subspaces it lives in. Two different kets must have the same system dimensions in order to be added together.\n\nIt is possible to normalize the ket vector after construction with the normalize! function.\n\nExample\n\njulia> ψ = normalize!(Ket([1,1]))\n2-d Schrodinger.Ket{Array{Float64,1},1} with space dimensions 2:\n0.71∠0°|0⟩ + 0.71∠0°|1⟩\n\n\n\n"
},

{
    "location": "api/quobj.html#Schrodinger.Operator",
    "page": "Quantum Object Types",
    "title": "Schrodinger.Operator",
    "category": "Type",
    "text": "Operator(B, dims=(size(B,1),))\n\nConstruct a linear operator from the matrix B. An N×N matrix will by default be assumed to describe an operator that acts on a single Hilbert space of dimension N. If the matrix represents a linear operator on a tensor product of Hilbert spaces, the dimensions can be defined manually by passing a tuple of subspace dimensions dims. In that case, prod(dims) must equal size(B,1).\n\nThe Operator type has two fields, data and dims, which store the matrix data and the subspace dimensions. An Operator, like a Ket or a Density matrix, is parameterized by the number of subspaces it lives in. Two different density matrices must have the same system dimensions in order to be added together. An Operator may or may not be Hermitian.\n\nExample\n\njulia> σ = Operator([0 -im ; im 0])\n2×2 Schrodinger.Operator{Array{Complex{Float64},2},1} with space dimensions 2:\n 0.0+0.0im  0.0-1.0im\n 0.0+1.0im  0.0+0.0im\n\n\n\n"
},

{
    "location": "api/quobj.html#Quantum-Object-Types-1",
    "page": "Quantum Object Types",
    "title": "Quantum Object Types",
    "category": "section",
    "text": "Modules = [Schrodinger]\nPages   = [\"quobj/types.jl\"]\nOrder   = [:type, :function]\nPrivate = false"
},

{
    "location": "api/states.html#",
    "page": "State Library",
    "title": "State Library",
    "category": "page",
    "text": "DocTestSetup  = quote\n    using Schrodinger\nend"
},

{
    "location": "api/states.html#Schrodinger.basis-Tuple{Integer,Integer}",
    "page": "State Library",
    "title": "Schrodinger.basis",
    "category": "Method",
    "text": "basis(N, n)\n\nGenerate a basis state (a.k.a. Fock or number state) ket n, in a Hilbert space of size N. Note that the size of the Hilbert space must be at least n+1. The function fock is an alias for basis.\n\nReturns a sparse vector.\n\nExample\n\njulia> ψ = basis(3,2)\n3-d Schrodinger.Ket{SparseVector{Float64,Int64},1} with space dimensions 3:\n1.00∠0°|2⟩\n\n\n\n"
},

{
    "location": "api/states.html#Schrodinger.coherent",
    "page": "State Library",
    "title": "Schrodinger.coherent",
    "category": "Function",
    "text": "coherent(N, α, analytic=false)\n\nGenerate a coherent state ket , in a Hilbert space of size N. To create a coherent density matrix, use the Density function: Density(coherent(N,n)).\n\nTwo methods can be used for generating a coherent state: via application of a displacment operator on a ground state (the default), or analytically, with the formula\n\n = e^-frac^22 sum_n=0^N-1 frac^nsqrtn n\n\nWhile the operator method will return a normalized ket, the analytic method will not. Both methods converge as N gets larger. The analytic method is also much faster, especially for large N.\n\nReturns a dense vector.\n\nExample\n\njulia> coherent(6,0.4+1im)\n6-d Schrodinger.Ket{Array{Complex{Float64},1},1} with space dimensions 6:\n0.60∠68°|1⟩ + 0.56∠0°|0⟩ + 0.46∠136°|2⟩ + 0.29∠-155°|3⟩ + 0.15∠-87°|4⟩\n\n\n\n"
},

{
    "location": "api/states.html#Schrodinger.maxmixed-Tuple{Integer}",
    "page": "State Library",
    "title": "Schrodinger.maxmixed",
    "category": "Method",
    "text": "maxmixed(N)\n\nGenerate a maximally mixed density matrix. The maximally mixed state is a mixture of basis states with uniform probability.\n\nReturns a sparse matrix.\n\nExample\n\njulia> maxmixed(4)\n4×4 Schrodinger.Density{SparseMatrixCSC{Float64,Int64},1} with space dimensions 4:\n 0.25  0.0   0.0   0.0\n 0.0   0.25  0.0   0.0\n 0.0   0.0   0.25  0.0\n 0.0   0.0   0.0   0.25\n\n\n\n"
},

{
    "location": "api/states.html#Schrodinger.thermal-Tuple{Integer,Real}",
    "page": "State Library",
    "title": "Schrodinger.thermal",
    "category": "Method",
    "text": "thermal(N, n)\n\nGenerate a thermal state density matrix _n with particle number n, in a Hilbert space of size N. A thermal state _n is a probabilistic mixture of basis states such that the expectation value of the number operator hatn is n. Note that this is true only if Nn. The returned density matrix is always normalized.\n\nReturns a sparse matrix.\n\nExample\n\njulia> N=5; n=0.2;\n\njulia> ρ = thermal(N,n)\n5×5 Schrodinger.Density{SparseMatrixCSC{Float64,Int64},1} with space dimensions 5:\n 0.833441  0.0       0.0        0.0         0.0\n 0.0       0.138907  0.0        0.0         0.0\n 0.0       0.0       0.0231511  0.0         0.0\n 0.0       0.0       0.0        0.00385852  0.0\n 0.0       0.0       0.0        0.0         0.000643087\n\njulia> expect(numberop(N),ρ)\n0.19935691318327978\n\n\n\n"
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
    "text": "DocTestSetup  = quote\n    using Schrodinger\nend"
},

{
    "location": "api/operators.html#Schrodinger.create-Tuple{Integer}",
    "page": "Operator Library",
    "title": "Schrodinger.create",
    "category": "Method",
    "text": "create(N)\n\nGenerate a quantum harmonic oscillator raising (creation) operator hata^ in a truncated Hilbert space of size N. Returns a sparse matrix.\n\nExample\n\njulia> create(4)\n4×4 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with space dimensions 4:\n 0.0  0.0      0.0      0.0\n 1.0  0.0      0.0      0.0\n 0.0  1.41421  0.0      0.0\n 0.0  0.0      1.73205  0.0\n\n\n\n"
},

{
    "location": "api/operators.html#Schrodinger.destroy-Tuple{Integer}",
    "page": "Operator Library",
    "title": "Schrodinger.destroy",
    "category": "Method",
    "text": "destroy(N)\n\nGenerate a quantum harmonic oscillator lowering (annihilation) operator hata in a truncated Hilbert space of size N. Returns a sparse matrix.\n\nExample\n\njulia> destroy(4)\n4×4 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with space dimensions 4:\n 0.0  1.0  0.0      0.0\n 0.0  0.0  1.41421  0.0\n 0.0  0.0  0.0      1.73205\n 0.0  0.0  0.0      0.0\n\n\n\n"
},

{
    "location": "api/operators.html#Schrodinger.displacementop-Tuple{Integer,Number}",
    "page": "Operator Library",
    "title": "Schrodinger.displacementop",
    "category": "Method",
    "text": "displacementop(N, α)\n\nGenerate a quantum harmonic oscillator displacement operator hatD() in a truncated Hilbert space of size N. Returns a dense matrix.\n\nhatD() = expleft(hata^ - ^*hataright)\n\nExample\n\njulia> displacementop(3,0.5im)\n3×3 Schrodinger.Operator{Array{Complex{Float64},2},1} with space dimensions 3:\n   0.88262+0.0im            0.0+0.439802im  -0.166001+0.0im\n       0.0+0.439802im  0.647859+0.0im             0.0+0.621974im\n -0.166001+0.0im            0.0+0.621974im    0.76524+0.0im\n\n\n\n"
},

{
    "location": "api/operators.html#Schrodinger.numberop-Tuple{Integer}",
    "page": "Operator Library",
    "title": "Schrodinger.numberop",
    "category": "Method",
    "text": "numberop(N)\n\nGenerate a number operator hatn in a Hilbert space of size N. Returns a sparse matrix.\n\nExample\n\njulia> numberop(4)\n4×4 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with space dimensions 4:\n 0.0  0.0  0.0  0.0\n 0.0  1.0  0.0  0.0\n 0.0  0.0  2.0  0.0\n 0.0  0.0  0.0  3.0\n\n\n\n"
},

{
    "location": "api/operators.html#Schrodinger.qeye",
    "page": "Operator Library",
    "title": "Schrodinger.qeye",
    "category": "Function",
    "text": "qeye(N, dims=(N,))\n\nGenerate an identity operator for a Hilbert space of size N. It is possible to specify the subspace dimensions with the dims argument. Returns a sparse matrix.\n\nExample\n\njulia> qeye(4,(2,2))\n4×4 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},2} with space dimensions 2⊗2:\n 1.0  0.0  0.0  0.0\n 0.0  1.0  0.0  0.0\n 0.0  0.0  1.0  0.0\n 0.0  0.0  0.0  1.0\n\n\n\n"
},

{
    "location": "api/operators.html#Schrodinger.qzero",
    "page": "Operator Library",
    "title": "Schrodinger.qzero",
    "category": "Function",
    "text": "qzero(N, dims=(N,))\n\nGenerate a zero operator for a Hilbert space of size N. It is possible to specify the subspace dimensions with the dims argument. Returns a sparse matrix.\n\nExample\n\njulia> qzero(4,(2,2))\n4×4 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},2} with space dimensions 2⊗2:\n 0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0\n\n\n\n"
},

{
    "location": "api/operators.html#Schrodinger.squeezeop-Tuple{Integer,Number}",
    "page": "Operator Library",
    "title": "Schrodinger.squeezeop",
    "category": "Method",
    "text": "squeezeop(N, z)\n\nGenerate a quantum harmonic oscillator squeeze operator hatS(z) in a truncated Hilbert space of size N. Returns a dense matrix.\n\nhatS(z) = expleft(frac12left(z^*hata^2 - zhata^2right)right)\n\nExample\n\njulia> squeezeop(3,0.5im)\n3×3 Schrodinger.Operator{Array{Complex{Float64},2},1} with space dimensions 3:\n 0.938148+0.0im       0.0+0.0im       0.0-0.346234im\n      0.0+0.0im       1.0+0.0im       0.0+0.0im\n      0.0-0.346234im  0.0+0.0im  0.938148+0.0im\n\n\n\n"
},

{
    "location": "api/operators.html#Operator-Library-1",
    "page": "Operator Library",
    "title": "Operator Library",
    "category": "section",
    "text": "Modules = [Schrodinger]\nPages   = [\"operators.jl\"]\nOrder   = [:function]\nPrivate = false"
},

{
    "location": "api/functions.html#",
    "page": "Function Library",
    "title": "Function Library",
    "category": "page",
    "text": "DocTestSetup  = quote\n    using Schrodinger\nend"
},

{
    "location": "api/functions.html#Schrodinger.ptrace-Tuple{Schrodinger.Density,Any}",
    "page": "Function Library",
    "title": "Schrodinger.ptrace",
    "category": "Method",
    "text": "ptrace(ρ, out)\n\nCompute the partial trace of a Density matrix or Operator ρ by tracing out the subsystems specified by out. Multiple subsystems can be traced out by passing a sorted tuple of subsystem indices.\n\nExample\n\nΦ₊ = normalize!(basis(2,0)⊗basis(2,0) + basis(2,1)⊗basis(2,1)) # Bell pair\nΨ₊ = normalize!(basis(2,0)⊗basis(2,1) + basis(2,1)⊗basis(2,0)) # Bell pair\nρ  = 0.25 * Density(Φ₊) + 0.75 * Density(Ψ₊) # density matrix\nptrace(ρ,2) # trace out qubit 2\n# output\n2×2 Schrodinger.Density{Array{Float64,2},1} with space dimensions 2:\n 0.5  0.0\n 0.0  0.5\n\n\n\n"
},

{
    "location": "api/functions.html#Schrodinger.ptrace-Tuple{Schrodinger.Ket,Any}",
    "page": "Function Library",
    "title": "Schrodinger.ptrace",
    "category": "Method",
    "text": "ptrace(ψ, out)\n\nCompute the partial trace of a state Ket or Bra ψ by tracing out the subsystems specified by out. Returns a density matrix. Multiple subsystems can be traced out by passing a sorted tuple of subsystem indices.\n\nExample\n\njulia> Φ₊ = normalize!(basis(2,0)⊗basis(2,0) + basis(2,1)⊗basis(2,1)) # Bell pair\n4-d Schrodinger.Ket{SparseVector{Float64,Int64},2} with space dimensions 2⊗2:\n0.71∠0°|0,0⟩ + 0.71∠0°|1,1⟩\n\njulia> ptrace(Φ₊,1) # trace out qubit 1\n2×2 Schrodinger.Density{Array{Float64,2},1} with space dimensions 2:\n 0.5  0.0\n 0.0  0.5\n\n\n\n"
},

{
    "location": "api/functions.html#Function-Library-1",
    "page": "Function Library",
    "title": "Function Library",
    "category": "section",
    "text": "Modules = [Schrodinger]\nPages   = [\"math/special.jl\",\"math/ptrace.jl\"]\nOrder   = [:function]\nPrivate = false"
},

]}
