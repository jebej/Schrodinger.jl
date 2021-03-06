# Julia Base math definitions translation
import Base: +, -, *, /, ^, real, imag, abs, abs2, round, sqrt, exp, log
import LinearAlgebra: kron, dot, transpose, adjoint, conj

# Additive identity and inverse
+(A::T) where {T<:QuObject} = A
-(A::T) where {T<:QuObject} = T(-A.data,A.dims)

# QuVector / Number algebra
+(x::Ket,b::Number) = Ket(x.data.+b,x.dims)
+(b::Number,x::Ket) = Ket(x.data.+b,x.dims)
-(x::Ket,b::Number) = Ket(x.data.-b,x.dims)
-(b::Number,x::Ket) = Ket(b.-x.data,x.dims)
*(x::Ket,b::Number) = Ket(x.data.*b,x.dims)
*(b::Number,x::Ket) = Ket(x.data.*b,x.dims)
/(x::Ket,b::Number) = Ket(x.data./b,x.dims)
+(x::Bra,b::Number) = Bra(x.data.+b,x.dims)
+(b::Number,x::Bra) = Bra(x.data.+b,x.dims)
-(x::Bra,b::Number) = Bra(x.data.-b,x.dims)
-(b::Number,x::Bra) = Bra(b.-x.data,x.dims)
*(x::Bra,b::Number) = Bra(x.data.*b,x.dims)
*(b::Number,x::Bra) = Bra(x.data.*b,x.dims)
/(x::Bra,b::Number) = Bra(x.data./b,x.dims)
^(x::T,b::Number) where {T<:QuVector} = throw(ArgumentError("cannot exponentiate $(tname(T))"))
^(x::T,b::Integer) where {T<:QuVector} = throw(ArgumentError("cannot exponentiate $(tname(T))"))

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

/(b::Number,x::T) where {T<:QuObject} = throw(ArgumentError("cannot divide number by $(tname(T))"))

# QuObject/QuObject Arithmetic
+(x::Ket,y::Ket) = (dimsmatch(x,y); Ket(x.data+y.data,x.dims))
-(x::Ket,y::Ket) = (dimsmatch(x,y); Ket(x.data-y.data,x.dims))
+(x::Bra,y::Bra) = (dimsmatch(x,y); Bra(x.data+y.data,x.dims))
-(x::Bra,y::Bra) = (dimsmatch(x,y); Bra(x.data-y.data,x.dims))
+(A::Operator,B::Operator) = (dimsmatch(A,B); Operator(A.data+B.data,A.dims))
-(A::Operator,B::Operator) = (dimsmatch(A,B); Operator(A.data-B.data,A.dims))
+(x::T1,y::T2) where {T1<:QuObject,T2<:QuObject} = throw(ArgumentError("cannot add $(tname(T1)) to $(tname(T2))"))
-(x::T1,y::T2) where {T1<:QuObject,T2<:QuObject} = throw(ArgumentError("cannot subtract $(tname(T2)) from $(tname(T1))"))

# QuObject/QuObject Algebra
*(ρ::Operator,σ::Operator) = (dimsmatch(ρ,σ); Operator(ρ.data*σ.data,ρ.dims))
*(σ::Operator,ψ::Ket) = (dimsmatch(σ,ψ); Ket(σ.data*ψ.data,ψ.dims))
*(ψ::Bra,ϕ::Ket) = (dimsmatch(ψ,ϕ); dotu(ψ.data,ϕ.data))
*(ψ::Bra,σ::Operator) = (dimsmatch(σ,ψ); Bra(transpose(σ.data)*ψ.data,ψ.dims))
*(ψ::Ket,ϕ::Bra) = (dimsmatch(ψ,ϕ); Operator(ψ.data*transpose(ϕ.data),ψ.dims))
dot(ρ::Operator,σ::Operator) = (dimsmatch(ρ,σ); dot(ρ.data,σ.data))
dot(ψ::Ket,ϕ::Ket) = (dimsmatch(ψ,ϕ); dot(ψ.data,ϕ.data))
# Rest is disallowed
*(x::T1,y::T2) where {T1<:QuObject,T2<:QuObject} = throw(ArgumentError("cannot multiply $(tname(T1)) with $(tname(T2))"))
/(x::T1,y::T2) where {T1<:QuObject,T2<:QuObject} = throw(ArgumentError("cannot divide $(tname(T1)) with $(tname(T2))"))

# Tensor Product
kron(x::Ket,y::Ket) = Ket(kron(x.data,y.data),(x.dims...,y.dims...))
kron(x::Bra,y::Bra) = Bra(kron(x.data,y.data),(x.dims...,y.dims...))
kron(A::Operator,B::Operator) = Operator(kron(A.data,B.data),(A.dims...,B.dims...))
kron(x::T1,y::T2) where {T1<:QuObject,T2<:QuObject} = throw(ArgumentError("cannot tensor $(tname(T1)) with $(tname(T2))"))

# Transposition and conjugation
adjoint(ψ::Ket) = Bra(ψ)
adjoint(ψ::Bra) = Ket(ψ)
adjoint(ρ::Operator) = Operator(adjoint(ρ.data),ρ.dims)
transpose(ρ::Operator) = Operator(transpose(ρ.data),ρ.dims)
conj(A::Ket) = Ket(conj.(A.data),A.dims)
conj(A::Bra) = Bra(conj.(A.data),A.dims)
conj(A::Operator) = Operator(conj.(A.data),A.dims)

# Math
sqrt(ρ::Operator) = Operator(sqrt(Matrix(ρ.data)),ρ.dims)
exp(ρ::Operator) = Operator(exp(Matrix(ρ.data)),ρ.dims)
log(ρ::Operator) = Operator(log(Matrix(ρ.data)),ρ.dims)
# TODO: julia is missing trig functions on matrices, we can do them via diagonalization

# Misc
round(x::T,args...;kwargs...) where {T<:QuObject} = T(round.(x.data,args...;kwargs...),x.dims)
real(x::QuObject) = real.(x.data)
imag(x::QuObject) = imag.(x.data)
abs(x::QuObject) = abs.(x.data)
abs2(x::QuObject) = abs2.(x.data)
