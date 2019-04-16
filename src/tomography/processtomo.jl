using Schrodinger, LinearAlgebra, BenchmarkTools

M = [ 47   736   395   358   421   383
     710   107   323   468   423   315
     338   428   128   646   429   338
     440   359   731    46   407   408
     403   340   349   400   753    36
     400   391   384   408    10   755]./(6*3*786) # normalize

function pdg_process_tomo(M,A,info=false)
    # projected gradient descent ML
    abs(sum(M)-1)<0.1 || throw(ArgumentError("measurement counts not normalized!"))
    d = 2
    C = Matrix{ComplexF64}(I,d^2,d^2)/d
    # objective and gradient functions
    f = C -> loglikelihood(M,C,A)
    ∇f = C -> loglikelihood_gradient(M,C,A)
    # metaparameters & initial variables
    μ = 3/2d^2; γ = 0.3
    c₁ = 1E6; c₂ = f(C)
    info && println("starting cost = $c₂")
    while c₁ - c₂ > 1E-10
        c₁, ∇c = c₂, ∇f(C)
        D = project_CPTP(C .- 1/μ .* ∇c) - C
        α = 1.0; Π = γ*real(D⋅∇c)
        while (c₂ = f(C.+α.*D)) > c₁ + α*Π
            α = α/2 # backtrack
        end
        @. C = C + α*D
    end
    info && println("final cost = $c₂")
    return C
end

apply_process(C::Operator,ψ::Ket) = apply_process(C,Operator(ψ))
apply_process(C::Operator,ρ::Operator) = ptrace((transpose(ρ)⊗qeye(size(ρ,1)))*C,1)

function gen_likelihood_model_1Q()
    # Generate the A matrix used to calculate likelihoods
    # The A matrix depends on the input states and measurement operators
    # We choose the 6 axial Bloch states as our inputs and measurements
    # |+⟩,|-⟩,|+i⟩,|-i⟩,|1⟩,|0⟩
    ρ = normalize!.(Operator.(Ket.([[1,1],[1,-1],[1,1im],[1,-1im],[0,1],[1,0]])))
    E = ρ ./ 3
    return mapreduce(transpose,vcat,[vec(full(transpose(ρ)⊗E)) for ρ ∈ ρ, E ∈ E])
end

function loglikelihood(M,C,A)
    # Binomial statistics for the measurement count probability, up to some irrelevant constant
    P = A*vec(C)
    return -real(transpose(vec(M))*log.(P))
end
function loglikelihood_gradient(M,C,A)
    P = A*vec(C)
    return unvec(-A'*(vec(M)./P))
end

unvec(C) = (d=isqrt(length(C)); reshape(C,(d,d)))

function project_CPTP(C)
    # generate helper objects
    MdagvecId,MdagM = TP_helper_matrices(C)
    D = Vector{real(eltype(C))}(undef,size(C,1))
    V = Matrix{eltype(C)}(undef,size(C))
    x₁ = copy(vec(C)); y₁ = zero(x₁);
    x₂ = copy(y₁); y₂ = copy(y₁)
    p = copy(y₁); q = copy(y₁)
    p_diff = 1.0; q_diff = 1.0
    # iterate through TP & CP projections
    while p_diff^2 + q_diff^2 + 2*abs(p⋅(x₂-x₁)) + 2*abs(q⋅(y₂-y₁)) > 1E-4
        y₂ = vec(project_TP(unvec(x₁+p),MdagvecId,MdagM))
        p_diff = norm(x₁-y₂,2)
        @. p = x₁ - y₂ + p
        x₂ = vec(project_CP(unvec(y₂+q),D,V))
        q_diff = norm(y₂-x₂,2)
        @. q = y₂ - x₂ + q
        x₁, x₂ = x₂, x₁
        y₁, y₂ = y₂, y₁
    end
    return unvec(x₁)
end

function project_CP(C,D,V)
    Schrodinger.hermfact!(D,V,Hermitian(C))
    D .= max.(D,0)
    return V*Diagonal(D)*V'
end


function project_TP(C,MdagvecId,MdagM)
    d⁻¹ = 1/isqrt(size(C,1))
    return vec(C) .- d⁻¹.*MdagM*vec(C) .+ d⁻¹.*MdagvecId
end

function TP_helper_matrices(C)
    d = isqrt(size(C,1))
    Id = Matrix{Int}(I,d,d); x = zeros(Int,1,d)
    # this can be done more efficiently, but prob doesn't matter
    M = sum((@inbounds x[i],x[mod1(i-1,d)] = 1,0; Id ⊗ x ⊗ Id ⊗ x) for i=1:d)
    return M'*vec(Id), M'M
end
