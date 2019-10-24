# Tomography Tests
using Schrodinger
using Schrodinger: mle_state_tomo, build_density_matrix, pdg_process_tomo,
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
ONEQUBIT = normalize!.([g,e1,(g+e1),(g-e1),(g+1im*e1),(g-1im*e1)])

@testset "State Tomography" begin
# Input states and POVM operators for state tomo
ψ = coherent(d,0.32)
E_m = Operator.(ONEQUBIT[[3,4,5,6,1,2]])
# Create the A matrix
A = @inferred state_likelihood_model(E_m./3)
# Simulate measurements
P = [real(expect(ψ,E)) for E ∈ E_m]
M = round.(Int,P*N .+ 25 .* randn(size(P)))./N
# Reconstruct state
ρ = Operator(build_density_matrix(mle_state_tomo(M,A).minimizer))
# Fidelity
F1 = @inferred fidelity2(ρ,ψ)
@test 1-F1 < 0.001
end

@testset "Process Tomography" begin
# Input states and POVM operators for process tomo
ρ_in = Operator.(ONEQUBIT[[3,4,5,6,1,2]])
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
    # Simulate tomography measurement by calculating theoretical probabilities
    # and drawing from a multinomial distribution
    #P = [real(expect(G*ρ*G',E)) for E ∈ E_m,ρ ∈ ρ_in]
    Ctrue = @inferred operator_to_choi(G) # ideal Choi matrix
    P = reshape(real(A*vec(data(Ctrue))),length(E_m),length(ρ_in))
    #M = mapslices(p->rand(Multinomial(N,p)),P,dims=2)./(d^2*N)
    M = .-(abs.(size(P,1)/2 .- abs.(P .+ randn(size(P))./2^12)) .- size(P,1)/2)./size(P,2)
    # Reconstruct the Choi matrix from the measurements
    Crec = Operator(@inferred(pdg_process_tomo(M,A,false)),(d,d))
    # Compare using the J distance, aka the trace norm of the difference
    @static if VERSION < v"1.0.0-"
        @test trace_norm(Ctrue-Crec)/2d < 0.004
    else
        @test @inferred(trace_norm(Ctrue-Crec))/2d < 0.004
    end
    # Compare with gate fidelity
    F1 = @inferred gate_fidelity_choi(Crec,G)
    F2 = mean(fidelity2(G*ψ,apply_process(Crec,ψ)) for ψ ∈ ONEQUBIT)
    @test 1-F1 < 0.002
    @test F1 ≈ F2 atol=1E-6

    # also test superoperator conversions (TODO: move somewhere else)
    @test Crec ≈ (Crec |> choi_to_kraus |> kraus_to_natural |> natural_to_choi)
    @test Crec ≈ (Crec |> choi_to_natural |> natural_to_kraus |> kraus_to_choi)
    @test Crec ≈ (Crec |> choi_to_chi |> chi_to_choi)
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
