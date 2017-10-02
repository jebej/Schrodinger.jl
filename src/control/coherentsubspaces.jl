immutable CoherentSubspaces{T} <: ObjectiveFunction
    δt::Float64 # timestep (fixed)
    Ut::Matrix{Complex128} # target unitary
    s::IntSet # coherent subspaces
    Hd::Matrix{T} # drift Hamiltonian
    Hc::Vector{Matrix{T}} # control Hamiltonian(s)
    u_last::Vector{Float64} # last control amplitudes, ordered by control, and then time
    # Propagator and eigendecompotion storage
    U::Vector{Matrix{Complex128}} # individual propagators
    X::Vector{Matrix{Complex128}} # forward cumulative propagators
    P::Vector{Matrix{Complex128}} # backwards cumulative propagators
    D::Vector{Vector{Float64}} # eigenvalues associated to each propagator
    V::Vector{Matrix{T}} # eigenvector matrix associated to each propagator
    # Workspace variables
    H::Matrix{T} # temporary storage for full Hamiltonian
    A::Matrix{Complex128} # temporary
    Jkj::Matrix{Complex128} # temporary
    cisDj::Vector{Complex128} # temporary
end

const IntCol = Union{AbstractVector{Int},IntSet,Set{Int},NTuple{N,Int} where N}
Base.convert(::Type{IntSet},r::IntCol) = IntSet(r)
Base.IntSet(elems::Vararg{Int}) = IntSet(elems)

function CoherentSubspaces(Ut::Operator,s::IntCol,Hd::Operator,Hc::Vector{<:Operator},t::Real,n::Integer)
    # TODO: Prob should do some error checking
    N = prod(dims(Hd))
    m = length(Hc)
    # Make sure we pass dense operators
    Ut_d = full(Ut)
    Hd_d = full(Hd)
    Hc_d = full.(Hc)
    # Storage for last control ampitudes, NaN for first run
    u_last = fill(NaN64,m*n)
    # Generate cache for various objects
    H = promote_type(typeof(Hd_d),eltype(Hc_d))(N,N)
    A = Matrix{Complex128}(N,N)
    U = [Matrix{Complex128}(N,N) for i=1:n]
    X = [Matrix{Complex128}(N,N) for i=1:n]
    P = [Matrix{Complex128}(N,N) for i=1:n]
    # For the exact derivative we need to store eigenvectors and eigenvalues
    D = [Vector{Float64}(N) for i=1:n] # eigenvalues are always real
    V = [similar(H) for i=1:n]
    # More temporary storage
    Jkj = Matrix{Complex128}(N,N)
    cisDj = Vector{Complex128}(N)
    return CoherentSubspaces{eltype(H)}(t/n,Ut_d,s,Hd_d,Hc_d,u_last,U,X,P,D,V,H,A,Jkj,cisDj)
end

function (O::CoherentSubspaces)(u)
    # Calculate forward propagators
    calc_fprops!(O.U,O.X,O.D,O.V,u,O.δt,O.Hd,O.Hc,O.H,O.u_last)
    Uf = O.X[end]; N = size(O.Ut,1)
    # Calculate coherent subspaces norm infidelity:
    # fₑ = 1 - |⟨Ut,Uf⟩|cs/N
    return 1 - inner_cs(O.Ut,Uf,O.s)/N
end

function (O::CoherentSubspaces)(::Val{:gradient},fp,u)
    # Calculate forward and backward propagators
    calc_fprops!(O.U,O.X,O.D,O.V,u,O.δt,O.Hd,O.Hc,O.H,O.u_last)
    calc_bprops!(O.P,O.U,O.Ut)
    n = length(O.U); m = length(O.Hc); N = size(O.Ut,1)
    # Calculate exact derivative of coherent subspaces norm infidelity
    x,y = _inner_cs_1(O.X[end],O.Ut,O.s)
    for j = 1:n
        O.cisDj .= cis.(O.D[j])
        for k = 1:m
            Jmat!(O.Jkj,O.Hc[k],O.cisDj,O.D[j],O.V[j],O.δt,O.A)
            j==1 ? copy!(O.A,O.Jkj) : A_mul_B!(O.A,O.Jkj,O.X[j-1])
            fp[(j-1)*m+k] = -inner_cs_grad(O.P[j],O.A,x,y,O.s)/N
        end
    end
    return fp
end


function inner_cs{T,S}(A::AbstractMatrix{T},B::AbstractMatrix{S},s::IntSet)
    # calculate |trace(A'*B)|cs (coherent subspaces)
    # this inner product cares only about relative phase between the subspaces contained in s
    x,y = _inner_cs_1(A,B,s)
    return sqrt(_inner_cs_2(x,y))
end

function inner_cs_grad{T,S}(Pj::AbstractMatrix{T},JXj::AbstractMatrix{S},x,y,s::IntSet)
    # calculate the partial derivative of |trace(A'*B)|cs (coherent subspaces)
    w,z = _inner_cs_1(Pj,JXj,s)
    res = real(w*x)
    @inbounds for i = 1:length(z)
        res += real(z[i]*y[i])
        res += real(w*normalize(x))*abs(y[i]) + real(z[i]*normalize(y[i]))*abs(x)
        for j = 1:i-1
            res += real(z[i]*normalize(y[i]))*abs(y[j]) + real(z[j]*normalize(y[j]))*abs(y[i])
        end
    end
    return res/sqrt(_inner_cs_2(x,y))
end


Base.normalize(z::Complex) = cis(angle(z))
Base.normalize(z::Real) = one(z)

function _inner_cs_1{T,S}(A::AbstractMatrix{T},B::AbstractMatrix{S},s::IntSet)
    m, n = size(A)
    size(B) == (m,n) || throw(DimensionMismatch("matrices must have the same dimensions"))
    x = zero(promote_type(T,S))
    y = Vector{typeof(x)}(n-length(s))
    t = 1
    @inbounds for j = 1:n
        a = zero(x)
        for i = 1:m
            a += conj(A[i,j])*B[i,j]
        end
        if j in s
            x += a
        else
            y[t] = a
            t += 1
        end
    end
    return x,y
end

function _inner_cs_2(x,y)
    res = abs2(x)
    @inbounds for i = 1:length(y)
        res += abs2(y[i])
        res += 2*abs(y[i])*abs(x)
        for j = 1:i-1
            res += 2*abs(y[i])*abs(y[j])
        end
    end
    return res
end
