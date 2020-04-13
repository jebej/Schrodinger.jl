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
    # create function with autodiff
    t0 = [(k=(i-2d)÷d+1; k*(2d-k)+1 == i ? 1/√d : 0.0) for i=1:d^2]
    od = OnceDifferentiable(t -> loglikelihood_statetomo(t,M,A), t0; autodiff = :forward)
    return optimize(od, t0, ConjugateGradient(), Optim.Options())
end

function loglikelihood_statetomo(t,M,A)
    # Binomial statistics for the measurement count probability, up to some irrelevant constant
    X = t_expectation_values(t,A)
    return sum(x -> ((p,m)=x; -m*log(real(p))), zip(X,M))
end

function t_expectation_values(t::AbstractVector,A::Matrix)
    # assemble ρ from the t-vector
    ρ = build_density_matrix(t)
    # calculate the probabilities by multiplying the vectorized density
    # operator with the model matrix
    return A*vec(ρ)
end

function build_density_matrix(t::AbstractVector)
    # reconstruct the density matrix from the t-parameters
    T = build_T_matrix(t)
    α = 1/sum(abs2,T)
    # ρ = α * T'*T
    ρ = T'*T
    ρ .*= α
    return ρ
end

function build_T_matrix(t::AbstractVector{S}) where S<:Real
    # build the T matrix from the t-parameters
    d = isqrt(length(t))
    d^2 == length(t) || throw(ArgumentError("invalid t-vector length!"))
    T = zeros(Complex{S},d,d)
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
    return T
end
