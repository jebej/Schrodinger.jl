# Misc. Tests
using Schrodinger
using Schrodinger: ComplexF64
using Compat.Test, Compat.LinearAlgebra, Compat.SparseArrays
include("reffuns.jl")
println("Testing Misc. Functions...")

@testset "Misc. Tests" begin

@testset "Partial Trace" begin
cases = [
    (rand(4,4), (2,), (2,2)),
    (rand(6,6), (1,), (2,3)),
    (rand(24,24), (1,3), (2,4,3)),
    (rand(36,36), (1,), (2,9,2)),
    (rand(48,48), (1,2), (3,4,2,2)),
    (rand(50,50), (2,), (5,5,2)),
    (rand(ComplexF64,72,72), (1,2,4), (3,4,2,3)),]

for (A, out, sysdims) in cases
    @test ptrace(A,out,sysdims)[1] ≈ ptrace_ref(A,out,sysdims)[1]
    #a = median(@benchmark(ptrace($A,$out,$sysdims)))
    #b = median(@benchmark(ptrace_ref($A,$out,$sysdims)))
    #println("ptrace speedup: $(ratio(a,b))")
end
end

end
