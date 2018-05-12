using Base.LinAlg.RealHermSymComplexHerm
import Base: IntSet, reverse

convert(::Type{IntSet},r::IntCol) = IntSet(r)
IntSet(elems::Vararg{Int}) = IntSet(elems)

reverse(n::Number) = n

tname(T::Type) = T.name.name

function randomsmooth(n::Integer,m::Integer)
    R = Matrix{Float64}(n,m)
    for i = 1:m
        R[:,i] = randomsmooth(n)
    end
    return R
end

function randomsmooth(n::Integer)
    x = linspace(0,1,n).'
    N = rand(4:7) # number of sines to use
    f = 2.0 .+ 10.*rand(N) # frequencies
    A = normalize!(1./f.*(rand(N).-0.5),2) # amplitudes
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
        return sqrt(sqrt(2π*n)*(n/e)^n*(1.0+1/12n+1/288n^2))
    else
        throw(ArgumentError("n must be larger than or equal to 0"))
    end
end

function rand_unitary{T<:LinAlg.BlasComplex}(::Type{T},n::Integer)
    # Generate a Haar distributed random unitary matrix
    # ref https://arxiv.org/pdf/math-ph/0609050.pdf
    A  = randn(T,n,n)
    qr = qrfact!(A)
    U  = full(qr[:Q])
    @inbounds for i = 1:n
        U[:,i] .*= sign(qr[:R][i,i]) # U = Q*Diagonal(diag(R)./abs.(diag(R)))
    end
    return U
end

Base.normalize(z::Complex) = z == 0 ? one(z) : z/abs(z)
Base.normalize(z::Real) = one(z)

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

expim(H::RealHermSymComplexHerm) = expim!(Matrix{Complex128}(H),copy(H))
expim!(H::RealHermSymComplexHerm) = expim!(Matrix{Complex128}(H),H)

function expim!(R::Matrix{Complex128},H::RealHermSymComplexHerm,Λ=Vector{Float64}(size(H,1)),U=Matrix(H),B=similar(R))
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
    A_mul_Bc!(R,B,U)
    return R
end

# index in a tensored system (warning, this function uses indices that start at 0)
function tindex{N}(inds,sysdims::NTuple{N,Int})
    i = 1 + inds[N]
    d = 1
    @inbounds for n = N-1:-1:1
        d *= sysdims[n+1]
        i += d * inds[n]
    end
    return i
end

# (warning, this function uses indices that start at 1)
function tindexr{N}(rinds,rsysdims::NTuple{N,Int})
    i = rinds[1]
    d = 1
    @inbounds for n = 2:N
        d *= rsysdims[n-1]
        i += d * (rinds[n]-1)
    end
    return i
end

revtuple{N}(t::NTuple{N,Any}) = ntuple(i->t[N+1-i],Val{N})
revinds{N}(t::NTuple{N,Any},ns::Int) = ntuple(i->ns+1-t[N+1-i],Val{N})
gettuple{N}(t1::NTuple,t2::NTuple{N,Any}) = ntuple(i->t1[t2[i]],Val{N})

ntuple_sans_m{n}(m,::Type{Val{n}}) = sorted_setdiff(ntuple(identity,Val{n}),(m,))

# Thanks to @mbauman on Discourse for the following
# https://discourse.julialang.org/t/type-stable-difference-of-tuples/3933/4
@inline function sorted_setdiff(t1::Tuple, t2::Tuple)
    if t1[1] == t2[1]
        sorted_setdiff(tail(t1), tail(t2))
    else
        (t1[1], sorted_setdiff(tail(t1), t2)...)
    end
end
@noinline sorted_setdiff(t1::Tuple{}, t2::Tuple) = throw(ArgumentError("duplicate or missing index $(t2[1])"))
sorted_setdiff(t1::Tuple, ::Tuple{}) = t1
sorted_setdiff(::Tuple{}, ::Tuple{}) = ()
