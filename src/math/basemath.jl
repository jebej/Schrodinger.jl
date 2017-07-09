# Julia Base math definitions translation
import Base: +, -, *, /, ^,
    A_mul_Bc, Ac_mul_B, Ac_mul_Bc, A_mul_Bt, kron, dot,
    ctranspose, transpose, conj, BLAS.dotu,
    sqrt, exp, log, real, imag, abs, abs2

# Additive identity and inverse
+{T<:QuObject}(A::T) = A
-{T<:QuObject}(A::T) = T(-A.data,A.dims)

# QuVector / Number algebra
+(x::Ket,b::Number) = Ket(x.data+b,x.dims)
+(b::Number,x::Ket) = Ket(x.data+b,x.dims)
-(x::Ket,b::Number) = Ket(x.data-b,x.dims)
-(b::Number,x::Ket) = Ket(b-x.data,x.dims)
*(x::Ket,b::Number) = Ket(x.data*b,x.dims)
*(b::Number,x::Ket) = Ket(x.data*b,x.dims)
/(x::Ket,b::Number) = Ket(x.data/b,x.dims)
/(b::Number,x::Ket) = throw(ArgumentError("cannot divide a number by a Ket"))
+(x::Bra,b::Number) = Bra(x.data+b,x.dims)
+(b::Number,x::Bra) = Bra(x.data+b,x.dims)
-(x::Bra,b::Number) = Bra(x.data-b,x.dims)
-(b::Number,x::Bra) = Bra(b-x.data,x.dims)
*(x::Bra,b::Number) = Bra(x.data*b,x.dims)
*(b::Number,x::Bra) = Bra(x.data*b,x.dims)
/(x::Bra,b::Number) = Bra(x.data/b,x.dims)
/(b::Number,x::Bra) = throw(ArgumentError("cannot divide a number by a Bra"))

# Operator/Number algebra
+(A::Operator,b::Number) = Operator(A.data+b*I,A.dims)
+(b::Number,A::Operator) = Operator(A.data+b*I,A.dims)
-(A::Operator,b::Number) = Operator(A.data-b*I,A.dims)
-(b::Number,A::Operator) = Operator(b*I-A.data,A.dims)
*(A::Operator,b::Number) = Operator(A.data*b,A.dims)
*(b::Number,A::Operator) = Operator(A.data*b,A.dims)
/(A::Operator,b::Number) = Operator(A.data/b,A.dims)
/(b::Number,A::Operator) = throw(ArgumentError("cannot divide a number by an Operator"))
^(A::Operator,b::Number) = Operator(A.data^b,A.dims)
^(A::Operator,b::Integer) = Operator(A.data^b,A.dims) # This method is needed for some reason

# QuObject/QuObject Algebra
+(x,y::Ket) = (dimsmatch(x,y);Ket(x.data+y.data,x.dims))
-(x::Ket,y::Ket) = (dimsmatch(x,y);Ket(x.data-y.data,x.dims))
+(x::Bra,y::Bra) = (dimsmatch(x,y);Bra(x.data+y.data,x.dims))
-(x::Bra,y::Bra) = (dimsmatch(x,y);Bra(x.data-y.data,x.dims))
+(A::Operator,B::Operator) = (dimsmatch(A,B);Operator(A.data+B.data,A.dims))
-(A::Operator,B::Operator) = (dimsmatch(A,B);Operator(A.data-B.data,A.dims))
kron(x::Ket,y::Ket) = Ket(kron(x.data,y.data),(x.dims...,y.dims...))
kron(x::Bra,y::Bra) = Bra(kron(x.data,y.data),(x.dims...,y.dims...))
kron(A::Operator,B::Operator) = Operator(kron(A.data,B.data),(A.dims...,B.dims...))

# Operator/Operator
for f in (:*, :A_mul_Bc, :Ac_mul_B, :Ac_mul_Bc)
    @eval ($f)(ρ::Operator,σ::Operator) = (dimsmatch(ρ,σ);Operator(($f)(ρ.data,σ.data),ρ.dims))
end

# Operator/Ket
*(σ::Operator,ψ::Ket) = (dimsmatch(σ,ψ);Ket(*(σ.data,ψ.data),σ.dims))
Ac_mul_B(σ::Operator,ψ::Ket) = (dimsmatch(σ,ψ);Ket(Ac_mul_B(σ.data,ψ.data),σ.dims))

# Bra/Operator
# remember that bras are actually stored as regular vectors (not row vectors)
*(ψ::Bra,σ::Operator) = (dimsmatch(σ,ψ);Bra(At_mul_B(σ.data,ψ.data),σ.dims))
A_mul_Bc(ψ::Bra,σ::Operator) = (dimsmatch(σ,ψ);Bra(*(conj.(σ.data),ψ.data),σ.dims))

# Ket/Ket and Bra
*(ψ::Bra,ϕ::Ket) = (dimsmatch(ψ,ϕ);dotu(ψ.data,ϕ.data))
dot(ψ::Ket,ϕ::Ket) = (dimsmatch(ψ,ϕ);dot(ψ.data,ϕ.data))
Ac_mul_B(ψ::Ket,ϕ::Ket) = (dimsmatch(ψ,ϕ);dot(ψ.data,ϕ.data))
A_mul_Bc(ψ::Ket,ϕ::Ket) = (dimsmatch(ψ,ϕ);Operator(A_mul_Bc(ψ.data,ϕ.data),ψ.dims))
*(ψ::Ket,ϕ::Bra) = (dimsmatch(ψ,ϕ);Operator(A_mul_Bt(ψ.data,ϕ.data),ψ.dims))

# Transposition and conjugation
ctranspose(ψ::Ket) = Bra(ψ)
ctranspose(ψ::Bra) = Ket(ψ)
ctranspose(ρ::Operator) = Operator(ctranspose(ρ.data),ρ.dims)
transpose(ρ::Operator) = Operator(transpose(ρ.data),ρ.dims)
conj(ρ::Operator) = Operator(conj.(ρ.data),ρ.dims)

# Math
sqrt(ρ::Operator) = Operator(sqrtm(ρ.data),ρ.dims)
exp(ρ::Operator) = Operator(expm(ρ.data),ρ.dims)
log(ρ::Operator) = Operator(logm(ρ.data),ρ.dims)
# TODO: julia is missing trig functions on matrices, we can do them via diagonalization

# Misc
real(x::Ket) = Ket(real.(x.data),x.dims)
imag(x::Ket) = Ket(imag.(x.data),x.dims)
abs(x::Ket) = Ket(abs.(x.data),x.dims)
abs2(x::Ket) = Ket(abs2.(x.data),x.dims)
real(x::Bra) = Bra(real.(x.data),x.dims)
imag(x::Bra) = Bra(imag.(x.data),x.dims)
abs(x::Bra) = Bra(abs.(x.data),x.dims)
abs2(x::Bra) = Bra(abs2.(x.data),x.dims)
real(x::Operator) = Operator(real.(x.data),x.dims)
imag(x::Operator) = Operator(imag.(x.data),x.dims)
abs(x::Operator) = Operator(abs.(x.data),x.dims)
abs2(x::Operator) = Operator(abs2.(x.data),x.dims)
