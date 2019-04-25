using Compat.LinearAlgebra: copytri!

function state_likelihood_model(Eₘ_list)
    # Generate the A matrix used to calculate likelihoods
    # The A matrix depends on the measurement operators
    sum(abs,data(sum(Eₘ_list))-I)<1E-14 ||
        throw(ArgumentError("Eₘ operators do not form a valid POVM!"))
    return transpose(mapreduce(vec∘full,hcat,Eₘ_list)::Matrix{ComplexF64})
end

function mle_state_tomo(M,A)
    (nm,d²) = size(A); d = isqrt(d²)
    length(M) == nm || throw(ArgumentError("number of measurements does not match A-matrix size!"))
    d^2 == d² || throw(ArgumentError("invalid number of columns in A-matrix!"))
    f = let ρ = Matrix{ComplexF64}(undef,d,d), X = Vector{ComplexF64}(undef,nm)
        t -> loglikehood_binomial(t,M,A,ρ,X,zero(ρ))
    end
    return optimize(f,[fill(1/√d,d); zeros(d²-d)])
end

function loglikehood_gaussian(t,M,A,ρ,X,T)
    # Gaussian statistics for the measurement count probability
    t_expectation_values!(X,t,A,ρ,T)
    return sum((real(p)-m)^2/2p for (p,m) in zip(X,M))
end

function loglikehood_binomial(t,M,A,ρ,X,T)
    # Binomial statistics for the measurement count probability, up to some irrelevant constant
    t_expectation_values!(X,t,A,ρ,T)
    return sum(-m*log(real(p)) for (p,m) in zip(X,M))
end

function t_expectation_values!(X,t::AbstractVector,A::Matrix,ρ::Matrix,T::Matrix)
    # calculate ρ from the t-vector
    if length(t) == 4
        @inbounds t₁,t₂,t₃,t₄ = t
        build_density_matrix!(ρ,t₁,t₂,t₃,t₄)
    else
        build_density_matrix!(ρ,t,T)
    end
    # calculate the probabilities by multiplying the vectorized density
    # operator with the A matrix
    mul!(X,A,vec(ρ))
end

function build_density_matrix(t::AbstractVector{S}) where S<:Real
    d = isqrt(length(t))
    d^2 == length(t) || throw(ArgumentError("invalid t-vector length!"))
    ρ = Matrix{complex(float(S))}(undef,d,d)
    return build_density_matrix!(ρ,t,zero(ρ))
end

function build_density_matrix!(ρ::Matrix,t₁::Real,t₂::Real,t₃::Real,t₄::Real)
    # reconstruct the density matrix from the t-parameters, special 2-d case
    t₁²,t₂²,t₃²,t₄² = t₁^2,t₂^2,t₃^2,t₄^2
    N = t₁²+t₂²+t₃²+t₄²
    ρ[1,1] = (t₁²+t₃²+t₄²)/N;  ρ[1,2] = t₂*(t₃-1im*t₄)/N
    ρ[2,1] = t₂*(t₃+1im*t₄)/N; ρ[2,2] = t₂²/N
    return ρ
end

function build_density_matrix!(ρ::Matrix,t::AbstractVector,T::Matrix)
    # reconstruct the density matrix from the t-parameters, general case
    # T must have its upper triangle zeroed out
    size(ρ) == size(T) || throw(DimensionMismatch("ρ & T do not have the same dimensions!"))
    d = isqrt(length(t))
    (d,d) == size(T) || throw(DimensionMismatch("ρ & T incompatible with the t-vector!"))
    for i = 1:d
        @inbounds T[i,i] = t[i]
    end
    r,c = 2,1
    for i = d+1:2:d^2
        @inbounds T[r,c] = complex(t[i],t[i+1])
        r += 1
        if r > d
            c = c + 1
            r = c + 1
        end
    end
    α = 1/real(dot(T,T))
    # Ac_mul_B!(ρ,T,T)
    # ρ .= α .* ρ
    copytri!(BLAS.herk!('U', 'C', α, T, zero(α), ρ), 'U', true)
    return ρ
end

M = [103, 633, 621, 115, 377, 354]

# only valid for a certain set of measurement, and for a single qubit
# here, the inputs are count values for |+⟩,|-⟩,|+i⟩,|-i⟩,|0⟩,|1⟩
# measurements, e.g. M[1] = Nm*⟨+|ρ|+⟩
# a physical density matrix can be represented by 4 parameters; assuming
# some measurement statistics, we minimize the loglikehood function
# t_expectation_values(t::Array) = (@inbounds t₁,t₂,t₃,t₄=t; t_expectation_values(t₁,t₂,t₃,t₄))
#
# function t_expectation_values(t₁,t₂,t₃,t₄)
#     # Calculate the expectation values for |+⟩,|-⟩,|+i⟩,|-i⟩,|0⟩,|1⟩ with the
#     # t-parameterized density matrix ρ(t) = T†*T/trace(T†*T), where
#     # T = [  t₁    0
#     #      t₃+it₄  t₂]
#     t₁²,t₂²,t₃²,t₄² = t₁^2,t₂^2,t₃^2,t₄^2
#     N = t₁²+t₂²+t₃²+t₄²
#     # ⟨+|ρ|+⟩, ⟨-|ρ|-⟩
#     Xp, Xm = 1/2 + t₂*t₃/N, 1/2 - t₂*t₃/N
#     # ⟨+i|ρ|+i⟩, ⟨-i|ρ|-i⟩
#     Yp, Ym = 1/2 + t₂*t₄/N, 1/2 - t₂*t₄/N
#     # ⟨0|ρ|0⟩, ⟨1|ρ|1⟩
#     Zp, Zm = (t₁²+t₃²+t₄²)/N, t₂²/N
#     return Xp,Xm,Yp,Ym,Zp,Zm
# end
