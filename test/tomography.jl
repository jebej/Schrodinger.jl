# Tomography Tests
using Schrodinger
using Schrodinger: pdg_process_tomo, apply_process, operator_to_choi, trace_norm,
    gate_fidelity_choi
using Compat.Test, Compat.LinearAlgebra, Compat.Random
if VERSION > v"0.7.0-"
    using Compat.Statistics: mean
end
println("Testing Process Tomography...")

@testset "Process Tomography" begin
# Build a few different variable for testing
d = 2; N = 10^5
g = basis(2,0); e1 = basis(2,1)
ONEQUBIT = normalize!.([g,e1,(g+e1),(g-e1),(g+1im*e1),(g-1im*e1)])

# Input states and POVM operators
ρ_in = Operator.(ONEQUBIT[[3,5,1,2]])
E_m  = [ρ_in; ([qeye(d)].-ρ_in)]./d^2
@test sum(E_m) == qeye(d)

# Create the A matrix
A = mapreduce(transpose,vcat,[vec(full(ρ⊗transpose(E))) for ρ ∈ ρ_in, E ∈ E_m])

# Test reconstruction with a couple of notable gates
T = Operator([1 0; 0 cis(π/4)])
H = Operator([1 1; 1 -1]/√2)
for G ∈ [σ0, σx, σy, σz, √σx, T, H]
    # Simulate tomography measurement by calculating theoretical probabilities
    # and drawing from a multinomial distribution
    #P = [real(expect(G*ρ*G',E)) for ρ ∈ ρ_in, E ∈ E_m]
    Strue = @inferred operator_to_choi(G)
    P = reshape(real(A*vec(data(Strue))),d^2,2*d^2)
    #M = mapslices(p->rand(Multinomial(N,p)),P,dims=2)./(d^2*N)
    M = abs.(P .+ randn(size(P))./2^12)./(d^2)
    # Reconstruct the Choi matrix from the measurements
    Srec = Operator(@inferred(pdg_process_tomo(M,A,false)),(2,2))
    # Compare using the J distance, aka the trace norm of the difference
    @static if VERSION < v"1.0.0-"
        @test (trace_norm(Strue-Srec))/2d < 0.004
    else
        @test @inferred(trace_norm(Strue-Srec))/2d < 0.004
    end
    # Compare with gate fidelity
    F1 = @inferred gate_fidelity_choi(Srec,G)
    F2 = mean(fidelity2(G*ψ,apply_process(Srec,ψ)) for ψ ∈ ONEQUBIT)
    @test 1-F1 < 0.002
    @test F1 ≈ F2 atol=1E-6
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
