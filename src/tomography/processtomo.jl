function process_likelihood_model(ρ_list::Vector,Eₘ_list::Vector)
    # Generate the A matrix used to calculate likelihoods
    # The A matrix depends on the input states and measurement operators
    dimsmatch(ρ_list,Eₘ_list)
    sum(abs,data(sum(Eₘ_list))-I)<1E-15 ||
        throw(ArgumentError("Eₘ operators do not form a valid POVM!"))
    prep = (Eₘ,ρ) -> transpose(vec(Array(ρ⊗transpose(Eₘ))))
    return reduce(vcat,vec(prep.(Eₘ_list,permutedims(ρ_list))))
end

# below is the actual projected gradient descent algorithm from Knee, et al.
# Quantum process tomography via completely positive and trace-preserving
# projection. Phys. Rev. A 98, 062336 (2018).

function pgd_process_tomo(M::Matrix{T}, A::Matrix; tol=1E-10, cptp_tol=1E-8, info=false, cbfun=nothing) where T
    # Choi process matrix reconstruction by maximum likelihood projected gradient descent
    size(A,1)==length(M) || throw(ArgumentError("A matrix inconsistent with number of measurements!"))
    abs(sum(M)-1)<1/4 || throw(ArgumentError("measurement counts not normalized!"))
    # infer space dimensions from A matrix
    d = isqrt(isqrt(size(A,2)))
    # initial Choi matrix guess, the identity map
    C = Matrix{Complex{T}}(I/d,d^2,d^2)
    # objective and gradient functions setup
    f = C -> loglikelihood(M,C,A)
    ∇f = C -> loglikelihood_gradient(M,C,A)
    # metaparameters & initial cost calculation
    μ⁻¹ = 2d^2/T(3); γ = T(0.3)
    c₁ = T(1E6); c₂ = f(C)
    info && println("start cost = $c₂")
    # iterate through projected gradient descent steps, with backtracking
    MdagM, 𝕀 = TP_helpers(C)
    stp = 0
    while c₁ - c₂ > tol
        stp += 1
        c₁, ∇C = c₂, ∇f(C)
        D = project_CPTP(C .- μ⁻¹.*∇C, MdagM, 𝕀, cptp_tol) - C
        α = one(T); Π = γ*real(D⋅∇C)
        while (c₂ = f(C .+ α.*D)) > c₁ + α*Π
            α = α/2 # backtrack
        end
        C .= C .+ α.*D
        if cbfun !== nothing # run callback function; e.g. to calc fidelity at each step
            cbfun(stp,c₂,C)
        end
    end
    info && println("final cost = $c₂, Δ = $(c₁-c₂), number of steps: $stp")
    return C
end

function loglikelihood(M::Matrix{T}, C::Matrix{Complex{T}}, A::Matrix) where {T<:Real}
    # Binomial statistics for the measurement count probability, up to some irrelevant constant
    logP = log.(max.(real.(A*vec(C)), T(1E-12)))
    return -(M⋅logP)
end

function loglikelihood_gradient(M::Matrix{T}, C::Matrix{Complex{T}}, A::Matrix) where {T<:Real}
    mMP = .-vec(M)./max.(real.(A*vec(C)), T(1E-12))
    return unvec(A'*mMP)
end

function project_CPTP(C::Matrix{Complex{T}}, MdagM, 𝕀, tol::Real=1E-8) where {T<:Real}
    # generate storage objects
    X₁ = copy(C); Y₁ = zero(X₁);
    X₂ = copy(Y₁); Y₂ = copy(Y₁)
    P = copy(Y₁); Q = copy(Y₁)
    ΔP = one(T); ΔQ = one(T)
    # iterate through TP & CP projections
    while ΔP^2 + ΔQ^2 + 2*abs(P⋅X₂-P⋅X₁) + 2*abs(Q⋅Y₂-Q⋅Y₁) > tol
        project_TP!(Y₂, X₁+P, MdagM, 𝕀)
        ΔP = norm2_diff(X₁,Y₂)
        P .= X₁ .- Y₂ .+ P
        project_CP!(X₂, Y₂+Q)
        ΔQ = norm2_diff(Y₂,X₂)
        Q .= Y₂ .- X₂ .+ Q
        X₁, X₂ = X₂, X₁
        Y₁, Y₂ = Y₂, Y₁
    end
    return X₁
end

function project_CP!(X::Matrix{Complex{T}}, C::Matrix{Complex{T}}) where {T<:Real}
    # Project the process onto the completely positive subspace by making the
    # Choi matrix positive semidefinite
    # We do this by taking the eigendecomposition, setting any negative
    # eigenvalues to 0, and reconstructing the Choi matrix
    D,V = hermfact!(Hermitian(hermitianize!(C)))
    D .= max.(D,zero(T))
    mul!(C,Diagonal(D),V')
    mul!(X,V,C)
end

function project_TP!(Y::Matrix{Complex{T}}, C::Matrix{Complex{T}}, MdagM::Matrix, 𝕀::Matrix) where {T<:Real}
    # Project the process onto the trace-preserving subspace
    d⁻¹ = one(T)/isqrt(size(C,1))
    mul!(vec(Y), MdagM, vec(C))
    Y .= C .- d⁻¹.*(Y.-𝕀)
end

function TP_helpers(C::Matrix{Complex{T}}) where {T<:Real}
    d = isqrt(size(C,1))
    𝕀 = Matrix{T}(I,d,d); k = zeros(Complex{T},1,d)
    # this could be done more efficiently, but it doesn't matter
    M = sum(i->(k[i]=1; k[mod1(i-1,d)]=0; (𝕀 ⊗ k) ⊗ (𝕀 ⊗ k)), 1:d)
    return M'*M, 𝕀⊗𝕀
end
