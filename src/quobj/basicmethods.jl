import Base: length, size, eltype, getindex, setindex!, similar, copy, hash,
    isequal, ==, convert, promote_rule, isapprox, show, @propagate_inbounds
import Compat.LinearAlgebra: checksquare, diag, complex, rank, normalize!,
    normalize, ishermitian, issymmetric, isdiag, triu, tril
import Compat: norm, opnorm
if VERSION < v"1.0.0-"
    import Base: full
end
# Special QuObject methods
data(A::QuObject) = A.data
dims(A::QuObject) = A.dims
isnormalized(A::QuVector) = norm(A) ≈ 1
isnormalized(A::QuMatrix) = trace(A) ≈ 1
isdensityop(A::QuMatrix) = isapproxhermitian(A) && isnormalized(A)
function checkdensityop(A::QuMatrix)
    herm = isapproxhermitian(A); norm = isnormalized(A)
    if !herm || !norm
        a = !herm ? "is not Hermitian" : ""
        b = !herm && !norm ? " and " : ""
        c = !norm ? "does not have unity trace (trace: $(real(trace(A))))" : ""
        throw(ArgumentError(string("not a density operator, matrix ",a,b,c)))
    end
end
dense(x::Ket) = Ket(Vector(x.data),x.dims)
dense(x::Bra) = Bra(Vector(x.data),x.dims)
dense(A::Operator) = Operator(Matrix(A.data),A.dims)
dimsmatch(A::QuObject,B::QuObject) =
    A.dims==B.dims || throw(DimensionMismatch("subspace dimensions do not match"))
function dimsmatch(As::AbstractVecOrTuple{QuObject},Bs::AbstractVecOrTuple{QuObject})
    for A ∈ As, B ∈ Bs; dimsmatch(A,B); end
end
dimsmatch(As::AbstractVecOrTuple{QuObject},B::QuObject) = dimsmatch(As,(B,))
dimsmatch(A::QuObject,Bs::AbstractVecOrTuple{QuObject}) = dimsmatch((A,),Bs)

# Tensored indexing methods
@propagate_inbounds getindex(A::QuVector,t::Tuple) =
    getindex(A.data,tensored_sub2ind(dims(A),t))
@propagate_inbounds getindex(A::QuMatrix,ti::Tuple,tj::Tuple) =
    getindex(A.data,tensored_sub2ind(dims(A),ti),tensored_sub2ind(dims(A),tj))
@propagate_inbounds setindex!(A::QuVector,v,t::Tuple) =
    setindex!(A.data,v,tensored_sub2ind(dims(A),t))
@propagate_inbounds setindex!(A::QuMatrix,v,ti::Tuple,tj::Tuple) =
    setindex!(A.data,v,tensored_sub2ind(dims(A),ti),tensored_sub2ind(dims(A),tj))

# Translate basic Base array methods to QuObjects
@propagate_inbounds getindex(A::QuObject,idx...) = getindex(A.data,idx...)
@propagate_inbounds setindex!(A::QuObject,v,idx...) = setindex!(A.data,v,idx...)
length(A::QuObject) = length(A.data)
size(A::QuObject,d) = size(A.data,d)
size(A::QuObject) = size(A.data)
checksquare(A::QuMatrix) = checksquare(A.data)
eltype(A::QuObject) = eltype(A.data)
diag(A::QuMatrix,k::Int=0) = diag(A.data,k)
triu(A::QuMatrix) = Operator(triu(A.data),A.dims)
tril(A::QuMatrix) = Operator(tril(A.data),A.dims)
full(A::QuObject) = Array(A.data)
complex(x::Ket) = Ket(complex(x.data),x.dims)
complex(x::Bra) = Bra(complex(x.data),x.dims)
complex(A::Operator) = Operator(complex(A.data),A.dims)
norm(x::QuVector,p::Integer=2) = norm(x.data,p)
norm(A::QuMatrix,p::Integer=2) = norm(A.data,p)
opnorm(A::QuMatrix,p::Integer=2) = opnorm(A.data,p)
trace(A::QuMatrix) = trace(A.data)
rank(A::QuMatrix) = rank(A.data)
normalize!(x::QuVector) = (normalize!(x.data,2);x)
normalize!(A::QuMatrix) = (rmul!(A.data,1/trace(A.data));A)
normalize(x::QuObject) = normalize!(copy(x))
scale!(A::QuObject,b::Number) = (rmul!(A.data,b);A)
scale(A::QuObject,b::Number) = scale!(copy(A))
ishermitian(A::Operator) = ishermitian(A.data)
isapproxhermitian(A::Operator) = isapproxhermitian(A.data)
issymmetric(A::QuMatrix) = issymmetric(A.data)
isunitary(A::QuMatrix) = isunitary(A.data)
isapproxunitary(A::QuMatrix) = isapproxunitary(A.data)
isdiag(A::QuMatrix) = isdiag(A.data)
similar(A::T) where {T<:QuObject}= T(similar(A.data),A.dims)
copy(A::T) where {T<:QuObject} = T(copy(A.data),A.dims)
hash(A::T,h::UInt) where {T<:QuObject} = hash(hash(A.data,hash(A.dims,hash(T))),h)
isequal(A::T,B::T) where {T<:QuObject} = isequal(A.dims,B.dims)&&isequal(A.data,B.data)
==(A::Ket,B::Ket) = isequal(A.dims,B.dims)&&(A.data==B.data)
==(A::Bra,B::Bra) = isequal(A.dims,B.dims)&&(A.data==B.data)
==(A::Operator,B::Operator) = isequal(A.dims,B.dims)&&(A.data==B.data)
isapprox(A::QuObject,B::QuObject;kwargs...) = isequal(A.dims,B.dims)&&isapprox(A.data,B.data;kwargs...)

# Conversion and promotion rules
convert(::Type{Ket{T,D}},x::Ket{S,D}) where {T,S,D} = Ket(convert(T,x.data),x.dims)
convert(::Type{Bra{T,D}},x::Bra{S,D}) where {T,S,D} = Bra(convert(T,x.data),x.dims)
convert(::Type{Operator{T,D}},A::Operator{S,D}) where {T,S,D} = Operator(convert(T,A.data),A.dims)
promote_rule(::Type{Ket{T,D}},::Type{Ket{S,D}}) where {T,S,D} = Ket{promote_type(T,S),D}
promote_rule(::Type{Bra{T,D}},::Type{Bra{S,D}}) where {T,S,D} = Bra{promote_type(T,S),D}
promote_rule(::Type{Operator{T,D}},::Type{Operator{S,D}}) where {T,S,D} = Operator{promote_type(T,S),D}

# Misc. conversions
import Compat.LinearAlgebra: Hermitian, Symmetric
Hermitian(A::Operator) = error("Use the `hermitian` function to make an Operator Hermitian")
hermitian(A::Operator) = Operator(Hermitian(A.data),A.dims)
Symmetric(A::Operator) = error("Use the `symmetric` function to make an Operator Symmetric")
symmetric(A::Operator) = Operator(Symmetric(A.data),A.dims)

# Show methods
function show(io::IO, A::T) where T<:QuMatrix
    n = size(A,1)
    print(io, "$n×$n $T with dimensions ", join(A.dims,'⊗'))
end
function show(io::IO, ::MIME"text/plain", A::T) where T<:QuMatrix
    show(io, A)
    println(io); print_array(IOContext(io,:compact=>true), A.data); println(io)
end
show(io::IO, ψ::QuVector) = print(io, braket(ψ))
function show(io::IO, ::MIME"text/plain", ψ::T) where T<:QuVector
    n = size(ψ,1)
    println(io, "$n-d $T with dimensions ", join(ψ.dims,'⊗'))
    println(io, braket(ψ))
end

# Bra-ket printing
function braket(ψ::QuVector, N::Int = 5)
    # N is the max number of kets/bras to show
    idx = findall(x->abs2(x)>2.5E-5,ψ.data) # keep amplitudes larger than 0.005
    isempty(idx) && return "0"
    perm = sortperm(ψ.data[idx], by=abs2, rev=true)
    idx = length(perm)>N ? idx[perm[1:N]] : idx[perm]
    coeffs = map(Array(ψ.data[idx])) do x
        @sprintf("%.2f∠%d°", abs(x), rad2deg(angle(x)))
    end
    labels = join.(tensored_ind2sub.((dims(ψ),), idx.-1), ",")
    return prettybraket(ψ,coeffs,labels) * (length(perm)>N ? " +…" : "")
end

prettybraket(::Ket,coeffs,labels) = join(string.(coeffs,"|",labels,"⟩")," + ")
prettybraket(::Bra,coeffs,labels) = join(string.(coeffs,"⟨",labels,"|")," + ")
