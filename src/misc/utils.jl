using LinearAlgebra: RealHermSymComplexHerm

tname(T::Type) = T.name.name

function randomsmooth(n::Integer,m::Integer)
    R = Matrix{Float64}(undef,n,m)
    for i = 1:m
        R[:,i] = randomsmooth(n)
    end
    return R
end

function randomsmooth(n::Integer)
    x = transpose(0:1/(n-1):1)
    N = rand(4:7) # number of sines to use
    f = 2.0 .+ 10 .* rand(N) # frequencies
    A = normalize!((rand(N).-0.5)./f,2) # amplitudes
    ϕ = 2π.*rand(N) # phase
    return vec(sum(A.*sin.(f.*x.+ϕ),1))
end

gaussian(x::Real,σ::Real) = exp(-0.5*(x/σ)^2)
gaussianprime(x::Real,σ::Real) = -exp(-0.5*(x/σ)^2)*x/σ^2

function gaussian_fwhm(x::Real,w)
    # w is the full width half maximum (FWHM)
    return exp(-4ln2*(x/w)^2)
end

function gaussianprime_fwhm(x::Real,w)
    # w is the full width half maximum (FWHM)
    return -x*ln2*exp2(3-4(x/w)^2)/w^2
end

function sqrtfact(n::Integer)
    if 0 <= n < 128
        return sqrtfact_table[n+1]
    elseif n >= 128 # Stirling series for values not stored in the table
        return sqrt(sqrt(2π*n)*(n/ℯ)^n*(1.0+1/12n+1/288n^2))
    else
        throw(ArgumentError("n must be larger than or equal to 0"))
    end
end

# pairwise sum of squared differences (SSD)
function norm2_diff(A::AbstractArray{T},B::AbstractArray{T}) where {T}
    n = length(A)
    n == length(B) || throw(DimensionMismatch())
    n == 0 ? abs2(zero(T)-zero(T)) : norm2_diff(A, B, 1, n)
end

function norm2_diff(A::AbstractArray, B::AbstractArray, i1::Integer, n::Integer)
    if n < 128
        @inbounds s = abs2(A[i1] - B[i1])
        for i in i1+1 : i1+n-1
            @inbounds s += abs2(A[i] - B[i])
        end
        return s
    else
        n2 = n÷2
        return norm2_diff(A, B, i1, n2) + norm2_diff(A, B, i1+n2, n-n2)
    end
end

function unwrap!(p)
    n,m = size(p,1),size(p,2)
    n < 2 && return p
    for j = 1:m, i = 2:n
        d = p[i,j] - p[i-1,j]
        if abs(d) > π
            p[i,j] -= floor((d+π)/2π) * 2π
        end
    end
    return p
end

expim(H::RealHermSymComplexHerm) = expim!(Matrix{complex(eltype(H))}(H),copy(H))
expim!(H::RealHermSymComplexHerm) = expim!(Matrix{complex(eltype(H))}(H),H)

function expim!(R::Matrix,H::RealHermSymComplexHerm,Λ=Vector{real(eltype(R))}(undef,size(H,1)),U=Matrix(H),B=similar(R))
    # First decompose H into U*Λ*Uᴴ
    #F = eigfact!(H); copy!(Λ,F.values); copy!(U,F.vectors)
    hermfact!(Λ,U,H)
    # Calculate the imaginary exponential of each eigenvalue and multiply
    # each column of U to obtain B = U*exp(iΛ)
    n = length(Λ)
    @inbounds for j = 1:n
        a = cis(Λ[j])
        @simd for i = 1:n
            B[i,j] = a * U[i,j]
        end
    end
    # Finally multiply B by Uᴴ to obtain U*exp(iΛ)*Uᴴ = exp(i*H)
    mul!(R,B,adjoint(U))
    return R
end

excited_sys(A::QuObject, subs::Vararg{Integer}) = excited_sys(dims(A), sys)
excited_sys(::Union{Dims{N},Val{N}}, subs::Vararg{Integer}) where N =
    ntuple(i->ifelse(i ∈ subs, 1, 0), Val(N))

lengthrange(a::Real,len::Integer) = UnitRange{typeof(a)}(a, oftype(a, a+len-1))

tensored_iterator(dims::Dims) = (revtuple(t) for t ∈ product(lengthrange.(0,revtuple(dims))...))
