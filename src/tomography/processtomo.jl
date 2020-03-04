function process_likelihood_model(ρ_list::Vector,Eₘ_list::Vector)
    # Generate the A matrix used to calculate likelihoods
    # The A matrix depends on the input states and measurement operators
    dimsmatch(ρ_list,Eₘ_list)
    sum(abs,data(sum(Eₘ_list))-I)<1E-15 ||
        throw(ArgumentError("Eₘ operators do not form a valid POVM!"))
    prep = (Eₘ,ρ) -> transpose(vec(full(ρ⊗transpose(Eₘ))))
    return reduce(vcat,vec(prep.(Eₘ_list,permutedims(ρ_list))))
end

# below is the actual projected gradient descent algorithm from Knee, et al.
# Quantum process tomography via completely positive and trace-preserving
# projection. Phys. Rev. A 98, 062336 (2018).

function pgd_process_tomo(M::Matrix, A::Matrix; tol=1E-10, cptp_tol=1E-8, info=false, cbfun=nothing)
    # Choi process matrix reconstruction by maximum likelihood projected gradient descent
    size(A,1)==length(M) || throw(ArgumentError("A matrix inconsistent with number of measurements!"))
    abs(sum(M)-1)<1/4 || throw(ArgumentError("measurement counts not normalized!"))
    # infer space dimensions from A matrix
    d = isqrt(isqrt(size(A,2)))
    # initial Choi matrix guess, the identity map
    C = Matrix{ComplexF64}(I/d,d^2,d^2)
    # objective and gradient functions setup
    f = C -> loglikelihood(M,C,A)
    ∇f = C -> loglikelihood_gradient(M,C,A)
    # metaparameters & initial cost calculation
    μ = 3/2d^2; γ = 0.3
    c₁ = 1E6; c₂ = f(C)
    info && println("start cost = $c₂")
    # iterate through projected gradient descent steps, with backtracking
    h = CPTP_helpers(C)
    stp = 0
    while c₁ - c₂ > tol
        stp += 1
        c₁, ∇C = c₂, ∇f(C)
        D = project_CPTP(C .- 1/μ.*∇C, h, cptp_tol) - C
        α = 1.0; Π = γ*real(D⋅∇C)
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

function loglikelihood(M::Matrix{T}, C::Matrix, A::Matrix) where {T<:Real}
    # Binomial statistics for the measurement count probability, up to some irrelevant constant
    logP = log.(max.(real.(A*vec(C)), 1E-16))
    return -real(vec(M)⋅logP)
end

function loglikelihood_gradient(M::Matrix, C::Matrix, A::Matrix)
    P = max.(real.(A*vec(C)), 1E-16)
    return unvec(-A'*(vec(M)./P))
end

function project_CPTP(C::Matrix, h::Tuple, tol::Real=1E-8)
    # generate storage objects
    X₁ = copy(vec(C)); Y₁ = zero(X₁);
    X₂ = copy(Y₁); Y₂ = copy(Y₁)
    P = copy(Y₁); Q = copy(Y₁)
    ΔP = 1.0; ΔQ = 1.0
    D, V, Mdagvec𝕀, MdagM = h
    # iterate through TP & CP projections
    while ΔP^2 + ΔQ^2 + 2*abs(P⋅X₂-P⋅X₁) + 2*abs(Q⋅Y₂-Q⋅Y₁) > tol
        project_TP!(Y₂, X₁+P, Mdagvec𝕀, MdagM)
        ΔP = norm2_diff(X₁,Y₂)
        P .= X₁ .- Y₂ .+ P
        project_CP!(X₂, Y₂+Q, D, V)
        ΔQ = norm2_diff(Y₂,X₂)
        Q .= Y₂ .- X₂ .+ Q
        X₁, X₂ = X₂, X₁
        Y₁, Y₂ = Y₂, Y₁
    end
    return unvec(X₁)
end

function project_CP!(X::Vector, vecC::Vector, D::Vector, V::Matrix)
    # Project the process onto the completely positive subspace by making the
    # Choi matrix positive semidefinite
    # We do this by taking the eigendecomposition, setting any negative
    # eigenvalues to 0, and reconstructing the Choi matrix
    C = unvec(vecC)
    hermitianize!(C)
    hermfact!(D,V,Hermitian(C))
    D .= max.(D,0)
    mul!(C,Diagonal(D),V')
    mul!(unvec(X),V,C)
end

function project_TP!(Y::Vector, vecC::Vector, Mdagvec𝕀::Vector, MdagM::Matrix)
    # Project the process onto the trace-preserving subspace
    d⁻¹ = 1/isqrt(isqrt(length(vecC)))
    mul!(Y, MdagM, vecC)
    Y .= vecC .- d⁻¹.*Y .+ d⁻¹.*Mdagvec𝕀
end

function CPTP_helpers(C::Matrix)
    D = Vector{real(eltype(C))}(undef,size(C,1))
    V = Matrix{eltype(C)}(undef,size(C))
    Mdagvec𝕀, MdagM = TP_helper_matrices(C)
    return D, V, Mdagvec𝕀, MdagM
end

function TP_helper_matrices(C::Matrix)
    d = isqrt(size(C,1))
    𝕀 = Matrix(1.0I,d,d); k = zeros(1,d)
    # this can be done more efficiently, but prob doesn't matter
    M = sum(i->(k[i]=1; k[mod1(i-1,d)]=0; (𝕀 ⊗ k) ⊗ (𝕀 ⊗ k)), 1:d)
    return M'*vec(𝕀), M'*M
end
