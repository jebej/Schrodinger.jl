# Superoperator Tests
using Schrodinger
using Test, LinearAlgebra
using Schrodinger: apply_process, operator_to_choi, trace_norm,
    gate_fidelity_choi, kraus_to_natural, natural_to_kraus, natural_to_choi,
    choi_to_natural, kraus_to_choi, choi_to_kraus, choi_to_chi, chi_to_choi

println("Testing Superoperator Functions...")

@testset "Representation Changes" begin
d = 4
for i = 1:20
    U = rand_unitary(d)
    U2 = qeye(d) ⊗ U
    @test operator_to_choi(U) ≈ sum(U2*(ket((i,i),(d,d))*ket((j,j),(d,d))')*U2' for i=0:d-1,j=0:d-1)
end
end
