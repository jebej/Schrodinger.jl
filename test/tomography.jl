# Tomography Tests
using Schrodinger, Distributions
using Schrodinger: mle_state_tomo, build_density_matrix, pgd_process_tomo,
    apply_process, operator_to_choi, trace_norm, gate_fidelity_choi,
    state_likelihood_model, process_likelihood_model

using Schrodinger: kraus_to_natural, natural_to_kraus, natural_to_choi,
    choi_to_natural, kraus_to_choi, choi_to_kraus, choi_to_chi, chi_to_choi

using Compat.Test, Compat.LinearAlgebra, Compat.Random
if VERSION > v"0.7.0-"
    using Compat.Statistics: mean
end
println("Testing Tomography...")

@testset "Tomography" begin
# Build a few different variable for testing
d = 2; N = 10^5
g = basis(2,0); e1 = basis(2,1)
ONEQUBIT = normalize!.([(g+e1),(g-e1),(g+1im*e1),(g-1im*e1),g,e1])

@testset "State Tomography" begin
# Input states and POVM operators for state tomo
ψ = coherent(d,0.32)
E_m = Operator.(ONEQUBIT)./3
# Create the A matrix
A = @inferred state_likelihood_model(E_m)
# Simulate tomography measurement by drawing from a binomial distribution
P = [real(expect(ψ,E)) for E ∈ E_m]
M = [rand(Binomial(N,p))/N for p in P]
# Reconstruct state
ρ = Operator(build_density_matrix(mle_state_tomo(M,A).minimizer))
# Fidelity
F1 = @inferred fidelity2(ρ,ψ)
@test 1-F1 < 0.001
end

@testset "Process Tomography" begin
# Input states and POVM operators for process tomo
ρ_in = Operator.(ONEQUBIT)
E_m  = ρ_in./3
@test sum(E_m) ≈ qeye(d)

# Create the likelihood model matrix ("A-matrix")
A = @inferred process_likelihood_model(ρ_in,E_m)

# Test reconstruction with a couple of notable gates
X½ = Operator([1 -im; -im 1]/√2)
Y½ = Operator([1 -1; 1 1]/√2)
T = Operator([1 0; 0 cis(π/4)])
H = Operator([1 1; 1 -1]/√2)
for G ∈ [σ0, σx, σy, σz, X½, Y½, T, H]
    # Simulate tomography measurement by drawing from a binomial distribution
    P = [real(expect(G*ρ*G',E)) for E ∈ E_m, ρ ∈ ρ_in]
    M = [rand(Binomial(N,p))/N for p in P]./size(P,2)
    # Reconstruct the Choi matrix from the measurements
    Crec = Operator(@inferred(pgd_process_tomo(M,A,info=false)),(d,d))
    # Compare using the J distance, aka the trace norm of the difference
    Ctrue = @inferred operator_to_choi(G) # ideal Choi matrix
    @test @inferred(trace_norm(Ctrue-Crec))/2d < 0.003
    # Compare with gate fidelity
    F1 = @inferred gate_fidelity_choi(Crec,G)
    F2 = mean(ψ->fidelity2(G*ψ,apply_process(Crec,ψ)), ONEQUBIT)
    @test 1-F1 < 0.002
    @test F1 ≈ F2

    # also test superoperator conversions (TODO: move somewhere else)
    @test Crec ≈ (Crec |> choi_to_kraus |> kraus_to_natural |> natural_to_choi)
    @test Crec ≈ (Crec |> choi_to_natural |> natural_to_kraus |> kraus_to_choi)
    @test Crec ≈ (Crec |> choi_to_chi |> chi_to_choi)
end

# Also test with some 2Q gates
d = 4
ρ_in = Operator.(vec(ONEQUBIT .⊗ permutedims(ONEQUBIT)))
E_m = ρ_in/9 # POVM
@test sum(E_m) ≈ qeye((2,2))
A = process_likelihood_model(ρ_in,E_m)
include(joinpath(@__DIR__,"twoqubit.jl"))

CNOT = Operator(
  [1 0 0 0
   0 1 0 0
   0 0 0 1
   0 0 1 0], (2,2))
CZ = Operator(
  [1 0 0 0
   0 1 0 0
   0 0 1 0
   0 0 0 -1], (2,2))
SWAP = Operator(
  [1 0 0 0
   0 0 1 0
   0 1 0 0
   0 0 0 1], (2,2))
iSWAP = Operator(
  [1    0  0    0
   0    0  1im  0
   0  1im  0    0
   0    0  0    1], (2,2))
sqrtSWAP = Operator(
  [1         0         0  0
   0  (1+im)/2  (1-im)/2  0
   0  (1-im)/2  (1+im)/2  0
   0         0         0  1], (2,2))
sqrtiSWAP = Operator(
  [1        0          0  0
   0    1/√(2)  1im/√(2)  0
   0  1im/√(2)    1/√(2)  0
   0         0         0  1], (2,2))
for G ∈ [CNOT, CZ, SWAP, sqrtSWAP, iSWAP, sqrtiSWAP]
   # Simulate tomography measurement by drawing from a binomial distribution
   P = [real(expect(G*ρ*G',E)) for E ∈ E_m,ρ ∈ ρ_in]
   M = [rand(Binomial(N,p))/N for p in P]./size(P,2)
   # Reconstruct the Choi matrix from the measurements
   Crec = Operator(@inferred(pgd_process_tomo(M,A,info=false)),(2,2,2,2))
   # Compare using the J distance, aka the trace norm of the difference
   Ctrue = @inferred operator_to_choi(G) # ideal Choi matrix
   @test @inferred(trace_norm(Ctrue-Crec))/2d < 0.003
   # Compare with gate fidelity
   F1 = @inferred gate_fidelity_choi(Crec,G)
   F2 = mean(ψ->fidelity2(G*ψ,apply_process(Crec,ψ)), TWOQUBIT)
   @test 1-F1 < 0.002
   @test F1 ≈ F2

   # also test superoperator conversions (TODO: move somewhere else)
   #@test Crec ≈ (Crec |> choi_to_kraus |> kraus_to_natural |> natural_to_choi)
   #@test Crec ≈ (Crec |> choi_to_natural |> natural_to_kraus |> kraus_to_choi)
   #@test Crec ≈ (Crec |> choi_to_chi |> chi_to_choi)
end

end


end

# the below is for historical testing purpose
# the measurements come from a poorly implemented Z gate
M = [ 47   736   395   358   421   383
     710   107   323   468   423   315
     338   428   128   646   429   338
     440   359   731    46   407   408
     403   340   349   400   753    36
     400   391   384   408    10   755]./(6*3*786)
# benchmark with CPTP tol of 1E-4
# BenchmarkTools.Trial:
#   memory estimate:  965.63 KiB
#   allocs estimate:  2349
#   --------------
#   minimum time:     873.043 μs (0.00% GC)
#   median time:      902.563 μs (0.00% GC)
#   mean time:        974.332 μs (6.05% GC)
#   maximum time:     2.970 ms (50.09% GC)
#   --------------
#   samples:          5116
#   evals/sample:     1
