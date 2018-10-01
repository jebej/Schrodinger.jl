using BenchmarkTools

H1 = rand(ComplexF64,4,4); H1 = Hermitian(H1+H1');
H2 = rand(Float64,4,4);    H2 = Hermitian(H2+H2');

H = H2

w, Z = eig(H)

w2, Z2 = hermfact!(similar(w), similar(Z), copy(H))

@assert w==w2 && Z==Z2

a = @benchmark eig(H) setup=(A=copy($H))
println("BASE eigendecomposition")
show(STDOUT,MIME"text/plain"(),a)
println("\n")
b = @benchmark hermfact!(w2, Z2, A) setup=(A=copy($H))
println("In-place eigendecomposition")
show(STDOUT,MIME"text/plain"(),b)
println()
