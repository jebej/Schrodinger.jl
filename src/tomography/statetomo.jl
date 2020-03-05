function state_likelihood_model(Eₘ_list)
    # Generate the model matrix used to calculate likelihoods
    # The model matrix depends on the measurement operators
    sum(abs,data(sum(Eₘ_list))-I)<1E-14 ||
        throw(ArgumentError("Eₘ operators do not form a valid POVM!"))
    return reduce(vcat,transpose.(vec.(full.(Eₘ_list))))
end

function mle_state_tomo(M,A)
    (nm,d²) = size(A); d = isqrt(d²)
    length(M) == nm || throw(ArgumentError("number of measurements does not match model matrix size!"))
    d^2 == d² || throw(ArgumentError("invalid number of columns in model matrix!"))
    f = let M = M, A = A, ρ = Matrix{ComplexF64}(undef,d,d), T = zero(ρ)
        t -> loglikehood_binomial(t,M,A,ρ,T)
    end
    return optimize(f, [(k=(i-2d)÷d+1; k*(2d-k)+1 == i ? 1/√d : 0.0) for i=1:d^2])
end

function loglikehood_gaussian(t,M,A,ρ,T)
    # Gaussian statistics for the measurement count probability
    X = t_expectation_values(t,A,ρ,T)
    return sum(x -> ((p,m)=x; (real(p)-m)^2/2p), zip(X,M))
end

function loglikehood_binomial(t,M,A,ρ,T)
    # Binomial statistics for the measurement count probability, up to some irrelevant constant
    X = t_expectation_values(t,A,ρ,T)
    return sum(x -> ((p,m)=x; -m*log(real(p))), zip(X,M))
end

function t_expectation_values(t::AbstractVector,A::Matrix,ρ::Matrix,T::Matrix)
    # assemble ρ from the t-vector
    build_density_matrix!(ρ,t,T)
    # calculate the probabilities by multiplying the vectorized density
    # operator with the model matrix
    return A*vec(ρ)
end

function build_density_matrix(t::AbstractVector{S}) where S<:Real
    d = isqrt(length(t))
    d^2 == length(t) || throw(ArgumentError("invalid t-vector length!"))
    ρ = Matrix{complex(float(S))}(undef,d,d)
    return build_density_matrix!(ρ,t,zero(ρ))
end

function build_density_matrix!(ρ::Matrix,t::AbstractVector,T::Matrix)
    # reconstruct the density matrix from the t-parameters, general case
    # IMPORTANT: T MUST have its upper triangle zeroed out
    d = isqrt(length(t))
    r,c,k = d,0,1
    @inbounds while k <= d^2
        r += 1
        if r > d
            c = c + 1
            r = c
            T[r,c] = complex(t[k])
            k += 1
        else
            T[r,c] = complex(t[k],t[k+1])
            k += 2
        end
    end
    # ρ = α * T'*T
    α = complex(1/sum(abs2,T))
    BLAS.gemm!('C', 'N', α, T, T, zero(α), ρ)
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
