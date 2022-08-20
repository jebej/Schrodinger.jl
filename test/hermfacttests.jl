using Schrodinger, LinearAlgebra, BenchmarkTools
using Schrodinger: hermfact!

d = 4

H1 = rand(ComplexF64,d,d); H1 = Hermitian(H1+H1');
H2 = rand(Float64,d,d);    H2 = Hermitian(H2+H2');

H = H1

w1, Z1 = eigen(H)

w2, Z2 = hermfact!(copy(H))


@assert w1==w2 && Z1==Z2


a = @benchmark eigen($H)
println("LinearAlgebra normal eigendecomposition:")
show(stdout, MIME"text/plain"(), a)
println("\n")

b = @benchmark eigen!(A) setup=(A=copy($H))
println("LinearAlgebra in-place eigendecomposition:")
show(stdout, MIME"text/plain"(), b)
println("\n")

c = @benchmark hermfact!(A) setup=(A=copy($H))
println("FastLapackInterface eigendecomposition:")
show(stdout, MIME"text/plain"(), c)
println("\n")
