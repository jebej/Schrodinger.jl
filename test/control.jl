using Schrodinger
using Test, LinearAlgebra, SparseArrays
using Schrodinger: gaussianpulse

println("Testing Optimal Control Functions...")

function opt_2Q_QFT(n)
    # Hamiltonian & controls
    Hd = 0.5 * (σx⊗σx + σy⊗σy + σz⊗σz)
    Hc = [0.5*σx⊗σ0, 0.5*σy⊗σ0, 0.5*σ0⊗σx, 0.5*σ0⊗σy]
    t = 6.0
    # Initial controls
    tvec = range(-t/2,t/2,length=n)
    ui = reduce(hcat, @. [2*gaussianpulse(tvec,[[t/4,t,2π/t,0,π]]),
                          2*gaussianpulse(tvec,[[t/4,t,2π/t,-π/2,π]]),
                          3*gaussianpulse(tvec,[[t/5,t,4π/t,0,π]]),
                          3*gaussianpulse(tvec,[[t/5,t,4π/t,-π/2,π]])])
    # Create objective function
    Ut = Operator([0.5  0.5    0.5  0.5
                   0.5  0.5im -0.5 -0.5im
                   0.5 -0.5    0.5 -0.5
                   0.5 -0.5im -0.5  0.5im],
                  (2,2))
    O = NormPSU(Ut,Hd,Hc,t,n)
    return grape(O,ui)
end

@testset "Optimal Control Tests" begin
    grape_res = opt_2Q_QFT(200)
    @test 1 - gate_fidelity(grape_res.Uf, grape_res.Ut) < 1E-8
end
