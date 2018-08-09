# QuObj Tests
using Base.Test, Schrodinger

@testset "Quantum Information" begin
# Build a few different variable for testing
N = 13
α = complex(1.1,-0.33)
g = basis(2,0)
e1 = basis(2,1)
numop = numberop(N)
const N! = normalize!

@testset "Expectation Value" begin
    # Test calculation of expectation valu
    @test expect(σz,g) == 1
    @test expect(e1,σz) == -1
    @test expect(σx,N!(g+e1)) ≈ 1
    @test expect(N!(g-e1),σx) ≈ -1
    @test expect(N!(g+1im*e1),σy) ≈ 1
    @test expect(σy,N!(g-1im*e1)) ≈ -1
    @test expect(coherent(N,α),destroy(N)) ≈ α
    @test expect(destroy(N),coherent(N,0.5α,true)) ≈ 0.5α
    @test [expect(basis(N,i),numop) for i=0:N-1] == collect(0:N-1)
    @test expect(N!(basis(N,2)+2*basis(N,3)),numop) == (2+2^2*3)/5
    @test expect(N!(Operator(basis(N,2))+2*Operator(basis(N,3))),numop) == (1*2+2*3)/3
    @test expect(numop,thermal(N,0.21)) ≈ 0.21
    @test expect(projectorop(N,[2,4]),maxmixed(N)) == 2/N
end

agf_test(U,V,S) = (X=U'V; sum(ψ->abs2(ψ⋅(X*ψ)),S)/length(S))
const ONEQUBIT = [g,e1,N!(g+e1),N!(g-e1),N!(g+1im*e1),N!(g-1im*e1)]
const TWOQUBIT = Ket.(N!.(Vector{Complex128}[
    # Separable states
    [0,0,1,0],[0,0,0,1],[0,0,1,1],[0,0,1,-1],[0,0,1,1im],[0,0,1,-1im],
    [1,0,-1,0],[0,1,0,-1],[1,1,-1,-1],[1,-1,-1,1],[1,1im,-1,-1im],[1,-1im,-1,1im],
    [1,0,-1im,0],[0,1,0,-1im],[1,1,-1im,-1im],[1,-1,-1im,1im],[1,1im,-1im,1],[1,-1im,-1im,-1],
    [1,0,1,0],[0,1,0,1],[1,1,1,1],[1,-1,1,-1],[1,1im,1,1im],[1,-1im,1,-1im],
    [1,0,1im,0],[0,1,0,1im],[1,1,1im,1im],[1,-1,1im,-1im],[1,1im,1im,-1],[1,-1im,1im,1],
    [1,0,0,0],[0,1,0,0],[1,1,0,0],[1,-1,0,0],[1,1im,0,0],[1,-1im,0,0],
    # Entangled states
    [0,1,1,0],[1,0,0,-1],[1,0,0,1],[0,1,-1,0],
    [1,0,0,1im],[0,1,1im,0],[0,1,-1im,0],[1,0,0,-1im],
    [1,1,1,-1],[1,1,-1,1],[1,-1,1,1],[1,-1,-1,-1],
    [1,1im,1,-1im],[1,1im,-1,1im],[1,-1im,1,1im],[1,-1im,-1,-1im],
    [1,1,1im,-1im],[1,1,-1im,1im],[1,-1,1im,1im],[1,-1,-1im,-1im],
    [1,1im,1im,1],[1,1im,-1im,-1],[1,-1im,1im,-1],[1,-1im,-1im,1],
    ]),((2,2),))

@testset "Fidelity" begin
    # Test calculation of fidelity
    @test fidelity(g,g) == 1
    @test fidelity2(g,g) == 1
    @test fidelity(g,e1) == 0
    @test fidelity2(g,e1) == 0
    @test fidelity(g,N!(g+e1)) == 1/√(2)
    @test fidelity2(g,N!(g+e1)) ≈ 1/2
    @test fidelity(Operator(coherent(N,α)),basis(N,3)) ≈ abs(coherent(N,α)[4])
    @test fidelity(basis(N,3),Operator(coherent(N,0.5α))) ≈ abs(coherent(N,0.5α)[4])
    @test fidelity2(Operator(coherent(N,0.5α)),basis(N,3)) ≈ abs2(coherent(N,0.5α)[4])
    @test fidelity2(basis(N,3),Operator(coherent(N,α))) ≈ abs2(coherent(N,α)[4])
    @test isapprox(fidelity(Operator(coherent(N,α)),Operator(basis(N,3))), abs(coherent(N,α)[4]), atol=1E-8)
    @test fidelity(Operator(basis(N,2)),Operator(coherent(N,1))) ≈ abs(coherent(N,1)[3])
    @test fidelity2(Operator(coherent(N,α,true)),Operator(basis(N,3))) ≈ abs2(coherent(N,α,true)[4])
    @test fidelity2(Operator(basis(N,2)),Operator(coherent(N,1))) ≈ abs2(coherent(N,1)[3])
    U1, V1 = rand_unitary(5), rand_unitary(5)
    @test entanglement_fidelity(U1,V1) ≈ (ε=qeye(5)⊗(U1'V1);ϕ=maxentangled(2,5);abs2(ϕ⋅(ε*ϕ)))
    U2, V2 = rand_unitary(2), rand_unitary(2)
    @test gate_fidelity(U2,V2) ≈ agf_test(U2,V2,ONEQUBIT)
    U3, V3 = rand_unitary(4,(2,2)), rand_unitary(4,(2,2))
    @test gate_fidelity(U3,V3) ≈ agf_test(U3,V3,TWOQUBIT)
end

@testset "Level Probabilities" begin
    # Test calculation of level probabilities
    @test levelprobs(g) == [1,0]
    @test levelprobs(e1) == [0,1]
    @test levelprobs(N!(e1+g)) ≈ [0.5,0.5]
    @test levelprobs(e1⊗g) == [0,0,1,0]
    @test levelprobs(N!(e1⊗g+g⊗e1)) ≈ [0,0.5,0.5,0]
    @test levelprobs(N!(e1⊗g+g⊗g),1:2) ≈ [[0.5,0.5],[1,0]]
end

function inner_cs_test(pe=0.1,ϕ=1.23)
    # σx on qubit 1 of 2-qubit system
    # no phase difference
    V1 = ((1-pe)*σx + pe*qeye(2)) .⊗ (Operator(basis(2,0)), Operator(basis(2,1)))
    # with phase difference
    V2 = (V1[1], cis(ϕ)*V1[2])
    # Ideal gate
    U = σx⊗Operator(basis(2,0)) + σx⊗Operator(basis(2,1))
    # inner products
    w = sqrt(sum(abs2,diag(sum(V2)'U)))/2
    x = abs(inner(sum(data.(V1)),data(U)))/4
    y = sum(abs,inner.(data.(V2),(data(U),)))/4
    z = abs(inner(sum(data.(V2)),data(U)))/4
    w,x,y,z
end


function anharm_NOT_RF(steps)
    dragpulse = Schrodinger.dragpulse_cos
    α = -0.4 * 2π # anharmonicity (GHz)
    Hd = α*Operator(basis(3,2)) # drift Hamiltonian
    σX, σY = 1/2*(create(3)+destroy(3)), 1im/2*(create(3)-destroy(3))
    λ = √2
    T = 5 # (ns)
    tspan = (-T/2,T/2)
    paramsx, paramsy = [0,T,0,α,λ,0,π], [0,T,0,α,λ,π/2,π]
    V = Operator(SchrodingerProp(Hd,((σX,dragpulse,paramsx),(σY,dragpulse,paramsy)),tspan,steps))
    #U = Operator([0 1 0; 1 0 0; 0 0 1]) # 3lvl NOT gate
    ## Fidelity
    #ϕ = N!(basis(3,0)⊗basis(3,0) + basis(3,1)⊗basis(3,1))
    #ϕ' * (qeye(3)⊗U'V) * Operator(ϕ) * (qeye(3)⊗V'U) * ϕ
    # Create objective function type
    #O = CoherentSubspaces(Ut,1:2,Hd,Hc,T,n) # care only about computational subspace
    #grape(O,ui)
    return V
end

end
