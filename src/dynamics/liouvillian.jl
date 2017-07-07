immutable Liouvillian{N,F,D} <: AbstractParameterizedFunction{true}
    L₀::SparseMatrixCSC{Complex128,Int}
    Lₙ::NTuple{N,SparseMatrixCSC{Complex128,Int}}
    fₙ::F
    pₙ::NTuple{N,Vector{Float64}}
    tmp::Vector{Complex128}
    dims::SDims{D}
    function (::Type{Liouvillian{N,D}}){N,D}(dims,L₀,Lₙ=(),fₙ=(),pₙ=())
        tmp = Vector{Complex128}(size(L₀,1))
        return new{N,typeof(fₙ),D}(L₀,Lₙ,fₙ,pₙ,tmp,dims)
    end
end

dimsmatch(L::Liouvillian,A::QuObject) = L.dims==dims(A) || throw(DimensionMismatch("subspace dimensions must match"))

function (L::Liouvillian{N,F,D}){N,F,D}(t,ψ,dψ) # Not sure why this function allocates 16 bytes...
    spmdv_mul!(dψ, L.L₀, ψ)
    for i = 1:N
        @inbounds applyfun!(dψ, L.Lₙ[i], L.fₙ[i], L.pₙ[i], t, ψ, L.tmp)
    end
end

function (L::Liouvillian{N,F,D}){N,F,D}(::Type{Val{:jac}},t,ψ,J)
    copy!(J, L.L₀)
    for i = 1:N
        @inbounds applyjac!(J, L.Lₙ[i], L.fₙ[i], L.pₙ[i])
    end
end

(p::Liouvillian)(::Type{Val{:hes}},t,ψ,H) = scale!(H, 0.0)

@inline function applyfun!(dψ,Lₙ,fₙ,pₙ,t,ψ,tmp)
    spmdv_mul!(tmp, Lₙ, ψ)
    scale!(tmp, fₙ(t,pₙ))
    dψ .+= tmp
end

@inline function applyjac!(J,Lₙ,fₙ,pₙ)
    J .+= scale!(Lₙ, fₙ(t,pₙ))
end

function spmdv_mul!(C::StridedVector, A::SparseMatrixCSC, B::StridedVector)
    nzv = A.nzval
    rv = A.rowval
    fill!(C, zero(eltype(C)))
    @inbounds for col = 1:A.n
        b = B[col]
        #if real(b) != 0.0 || imag(b) != 0.0
            for j = A.colptr[col]:(A.colptr[col+1]-1)
                C[rv[j]] += nzv[j]*b
            end
        #end
    end
    return C
end
