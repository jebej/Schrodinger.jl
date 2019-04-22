function gen_likelihood_model_1Q()
    # Generate the A matrix used to calculate likelihoods
    # The A matrix depends on the input states and measurement operators
    # We choose the 6 axial Bloch states as our inputs and measurements
    # |+⟩,|-⟩,|+i⟩,|-i⟩,|1⟩,|0⟩
    ρ = normalize!.(Operator.(Ket.([[1,1],[1,-1],[1,1im],[1,-1im],[1,0],[0,1]])))
    E = ρ ./ 3
    return mapreduce(transpose,vcat,[vec(full(transpose(ρ)⊗E)) for ρ ∈ ρ, E ∈ E])
end

# below is the actual projected gradient descent algorithm from Knee, et al.
# Quantum process tomography via completely positive and trace-preserving
# projection. Phys. Rev. A 98, 062336 (2018).

function pdg_process_tomo(M,A,info=false)
    # Choi process matrix reconstruction with maximum likelihood projected
    # gradient descent
    size(A,1)==length(M) || throw(ArgumentError("A matrix inconsistent with number of measurements!"))
    abs(sum(M)-1)<0.1 || throw(ArgumentError("measurement counts not normalized!"))
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
    while c₁ - c₂ > 1E-10
        c₁, ∇c = c₂, ∇f(C)
        D = project_CPTP(C.-1/μ.*∇c, h) - C
        α = 1.0; Π = γ*real(D⋅∇c)
        while (c₂ = f(C.+α.*D)) > c₁ + α*Π
            α = α/2 # backtrack
        end
        @. C = C + α*D
    end
    info && println("final cost = $c₂, Δc = $(c₁-c₂)")
    return C
end

function loglikelihood(M,C,A)
    # Binomial statistics for the measurement count probability, up to some irrelevant constant
    P = real.(A*vec(C))
    P .= max.(P,1E-16)
    return -real(transpose(vec(M))*log.(P))
end

function loglikelihood_gradient(M,C,A)
    P = real.(A*vec(C))
    P .= max.(P,1E-16)
    return unvec(-A'*(vec(M)./P))
end

function project_CPTP(C,h)
    # generate storage objects
    x₁ = copy(vec(C)); y₁ = zero(x₁);
    x₂ = copy(y₁); y₂ = copy(y₁)
    p = copy(y₁); q = copy(y₁)
    p_diff = 1.0; q_diff = 1.0
    D,V,Mdagvec𝕀,MdagM = h
    # iterate through TP & CP projections
    while p_diff^2 + q_diff^2 + 2*abs(p⋅(x₂-x₁)) + 2*abs(q⋅(y₂-y₁)) > 1E-6
        y₂ = project_TP(x₁+p,Mdagvec𝕀,MdagM)
        p_diff = norm(x₁-y₂,2)
        @. p = x₁ - y₂ + p
        x₂ = project_CP(y₂+q,D,V)
        q_diff = norm(y₂-x₂,2)
        @. q = y₂ - x₂ + q
        x₁, x₂ = x₂, x₁
        y₁, y₂ = y₂, y₁
    end
    return unvec(x₁)
end

function project_CP(vecC,D,V)
    # Project the process onto the completely positive subspace by making the
    # Choi matrix positive semidefinite
    # We do this by taking the eigendecomposition, setting any negative
    # eigenvalues to 0, and reconstructing the Choi matrix
    @static if VERSION < v"0.7.0-"
        H = unvec(vecC)
        @inbounds for i = 1:size(H,1); H[i,i] = real(H[i,i]); end
        hermfact!(D,V,Hermitian(H))
    else
        hermfact!(D,V,Hermitian(unvec(vecC)))
    end
    D .= max.(D,0)
    return vec(V*Diagonal(D)*V')
end

function project_TP(vecC,Mdagvec𝕀,MdagM)
    # Project the process onto the trace-preserving subspace
    d⁻¹ = 1/isqrt(isqrt(length(vecC)))
    return vecC .- d⁻¹.*MdagM*vecC .+ d⁻¹.*Mdagvec𝕀
end

function CPTP_helpers(C)
    D = Vector{real(eltype(C))}(undef,size(C,1))
    V = Matrix{eltype(C)}(undef,size(C))
    Mdagvec𝕀,MdagM = TP_helper_matrices(C)
    return D,V,Mdagvec𝕀,MdagM
end

function TP_helper_matrices(C)
    d = isqrt(size(C,1))
    𝕀 = Matrix{Int8}(I,d,d); k = zeros(Int8,1,d)
    # this can be done more efficiently, but prob doesn't matter
    M = sum((@inbounds k[i],k[mod1(i-1,d)] = 1,0; 𝕀 ⊗ k ⊗ 𝕀 ⊗ k) for i=1:d)
    return M'*vec(𝕀), M'*M
end
