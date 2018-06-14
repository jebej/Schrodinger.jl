# Julia Base math definitions translation
import Base: +, -, *, /, ^,
    A_mul_Bc, Ac_mul_B, Ac_mul_Bc, A_mul_Bt, kron, dot, vecdot,
    ctranspose, transpose, conj, BLAS.dotu,
    sqrtm, expm, logm, real, imag, abs, abs2, round

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
+(x::Bra,b::Number) = Bra(x.data+b,x.dims)
+(b::Number,x::Bra) = Bra(x.data+b,x.dims)
-(x::Bra,b::Number) = Bra(x.data-b,x.dims)
-(b::Number,x::Bra) = Bra(b-x.data,x.dims)
*(x::Bra,b::Number) = Bra(x.data*b,x.dims)
*(b::Number,x::Bra) = Bra(x.data*b,x.dims)
/(x::Bra,b::Number) = Bra(x.data/b,x.dims)
^{T<:QuVector}(x::T,b::Number) = throw(ArgumentError("cannot exponentiate $(tname(T))"))
^{T<:QuVector}(x::T,b::Integer) = throw(ArgumentError("cannot exponentiate $(tname(T))"))

# Operator/Number algebra
+(A::Operator,b::Number) = Operator(A.data+b*I,A.dims)
+(b::Number,A::Operator) = Operator(A.data+b*I,A.dims)
-(A::Operator,b::Number) = Operator(A.data-b*I,A.dims)
-(b::Number,A::Operator) = Operator(b*I-A.data,A.dims)
*(A::Operator,b::Number) = Operator(A.data*b,A.dims)
*(b::Number,A::Operator) = Operator(A.data*b,A.dims)
/(A::Operator,b::Number) = Operator(A.data/b,A.dims)
^(A::Operator,b::Number) = Operator(A.data^b,A.dims)
^(A::Operator,b::Integer) = Operator(A.data^b,A.dims) # This method is needed for some reason

/{T<:QuObject}(b::Number,x::T) = throw(ArgumentError("cannot divide number by $(tname(T))"))

# QuObject/QuObject Algebra
+(x::Ket,y::Ket) = (dimsmatch(x,y);Ket(x.data+y.data,x.dims))
-(x::Ket,y::Ket) = (dimsmatch(x,y);Ket(x.data-y.data,x.dims))
+(x::Bra,y::Bra) = (dimsmatch(x,y);Bra(x.data+y.data,x.dims))
-(x::Bra,y::Bra) = (dimsmatch(x,y);Bra(x.data-y.data,x.dims))
+(A::Operator,B::Operator) = (dimsmatch(A,B);Operator(A.data+B.data,A.dims))
-(A::Operator,B::Operator) = (dimsmatch(A,B);Operator(A.data-B.data,A.dims))
+{T1<:QuObject,T2<:QuObject}(x::T1,y::T2) = throw(ArgumentError("cannot add $(tname(T1)) to $(tname(T2))"))
-{T1<:QuObject,T2<:QuObject}(x::T1,y::T2) = throw(ArgumentError("cannot subtract $(tname(T2)) from $(tname(T1))"))

# Operator/Operator
*(ρ::Operator,σ::Operator) = (dimsmatch(ρ,σ);Operator(*(ρ.data,σ.data),ρ.dims))
A_mul_Bc(ρ::Operator,σ::Operator) = (dimsmatch(ρ,σ);Operator(A_mul_Bc(ρ.data,σ.data),ρ.dims))
Ac_mul_B(ρ::Operator,σ::Operator) = (dimsmatch(ρ,σ);Operator(Ac_mul_B(ρ.data,σ.data),ρ.dims))
Ac_mul_Bc(ρ::Operator,σ::Operator) = (dimsmatch(ρ,σ);Operator(Ac_mul_Bc(ρ.data,σ.data),ρ.dims))
vecdot(ρ::Operator,σ::Operator) = (dimsmatch(ρ,σ);vecdot(ρ.data,σ.data))

# Operator/Ket
*(σ::Operator,ψ::Ket) = (dimsmatch(σ,ψ);Ket(*(σ.data,ψ.data),ψ.dims))
Ac_mul_B(σ::Operator,ψ::Ket) = (dimsmatch(σ,ψ);Ket(Ac_mul_B(σ.data,ψ.data),ψ.dims))
Ac_mul_B(ψ::Ket,σ::Operator) = (dimsmatch(σ,ψ);Bra(vec(Ac_mul_B(ψ.data,σ.data)),ψ.dims))
Ac_mul_Bc(ψ::Ket,σ::Operator) = (dimsmatch(σ,ψ);Bra(vec(Ac_mul_Bc(ψ.data,σ.data)),ψ.dims))

# Bra/Operator
# remember that bras are actually stored as regular vectors (not row vectors)
*(ψ::Bra,σ::Operator) = (dimsmatch(σ,ψ);Bra(At_mul_B(σ.data,ψ.data),ψ.dims))
A_mul_Bc(ψ::Bra,σ::Operator) = (dimsmatch(σ,ψ);Bra(*(conj.(σ.data),ψ.data),ψ.dims))
A_mul_Bc(σ::Operator,ψ::Bra) = (dimsmatch(σ,ψ);Ket(*(σ.data,conj.(ψ.data)),ψ.dims))
Ac_mul_Bc(σ::Operator,ψ::Bra) = (dimsmatch(σ,ψ);Ket(Ac_mul_B(σ.data,conj.(ψ.data)),ψ.dims))

# Ket/Ket and Bra
dotu{T<:Real}(x::AbstractVector{T},y::AbstractVector{T}) = dot(x,y)
*(ψ::Bra,ϕ::Ket) = (dimsmatch(ψ,ϕ);dotu(ψ.data,ϕ.data))
dot(ψ::Ket,ϕ::Ket) = (dimsmatch(ψ,ϕ);dot(ψ.data,ϕ.data))
Ac_mul_B(ψ::Ket,ϕ::Ket) = dot(ψ,ϕ)
A_mul_Bc(ψ::Ket,ϕ::Ket) = (dimsmatch(ψ,ϕ);Operator(A_mul_Bc(ψ.data,ϕ.data),ψ.dims))
*(ψ::Ket,ϕ::Bra) = (dimsmatch(ψ,ϕ);Operator(A_mul_Bt(ψ.data,ϕ.data),ψ.dims))

# Rest is disallowed
*{T1<:QuObject,T2<:QuObject}(x::T1,y::T2) = throw(ArgumentError("cannot multiply $(tname(T1)) with $(tname(T2))"))
/{T1<:QuObject,T2<:QuObject}(x::T1,y::T2) = throw(ArgumentError("cannot divide $(tname(T1)) with $(tname(T2))"))

# Tensor Product
kron(x::Ket,y::Ket) = Ket(kron(x.data,y.data),(x.dims...,y.dims...))
kron(x::Bra,y::Bra) = Bra(kron(x.data,y.data),(x.dims...,y.dims...))
kron(A::Operator,B::Operator) = Operator(kron(A.data,B.data),(A.dims...,B.dims...))
kron{T1<:QuObject,T2<:QuObject}(x::T1,y::T2) = throw(ArgumentError("cannot tensor $(tname(T1)) with $(tname(T2))"))

# Transposition and conjugation
ctranspose(ψ::Ket) = Bra(ψ)
ctranspose(ψ::Bra) = Ket(ψ)
ctranspose(ρ::Operator) = Operator(ctranspose(ρ.data),ρ.dims)
transpose(ρ::Operator) = Operator(transpose(ρ.data),ρ.dims)
conj{T<:QuObject}(A::T) = T(conj.(A.data),A.dims)

# Math
sqrtm(ρ::Operator) = Operator(sqrtm(full(ρ)),ρ.dims)
expm(ρ::Operator) = Operator(expm(full(ρ)),ρ.dims)
logm(ρ::Operator) = Operator(logm(full(ρ)),ρ.dims)
# TODO: julia is missing trig functions on matrices, we can do them via diagonalization

# Misc
round(x::Operator,args...) = Operator(round.(x.data,args...),x.dims)
real(x::QuObject) = real.(x.data)
imag(x::QuObject) = imag.(x.data)
abs(x::QuObject) = abs.(x.data)
abs2(x::QuObject) = abs2.(x.data)
