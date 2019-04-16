
function super_basic_tomo(mx,my,mz)
    # here, the inputs are measurements in the Pauli basis
    (mx^2+my^2+mz^2)<=1 || throw(ArgumentError("invalid measurement values!"))
    ρ = [     1 + mz  mx - 1im*my
         mx + 1im*my       1 - mz]/2
    return Operator(ρ,(2,))
end

function basic_mle_tomo(M,Nm::Integer)
    # here, the inputs are count values for |+⟩,|-⟩,|+i⟩,|-i⟩,|1⟩,|0⟩
    # measurements, e.g. M[1] = Nm*⟨+|ρ|+⟩
    # a physical density matrix can be represented by 4 parameters; assuming
    # some measurement statistics, we minimize the loglikehood function
    length(M)==6 || throw(ArgumentError("measurement array must contain 6 values!"))
    f = T -> loglikehood_binomial(M,T,Nm)
    return optimize(f,[1/√2,1/√2,0,0])
end

function reconstruct_density(t₁,t₂,t₃,t₄)
    # reconstruct the density matrix from the t-parameters
    t₁²,t₂²,t₃²,t₄² = t₁^2,t₂^2,t₃^2,t₄^2
    N = t₁²+t₂²+t₃²+t₄²
    ρ = 1/N*[  t₁²+t₃²+t₄²   t₂*(t₃-1im*t₄)
             t₂*(t₃+1im*t₄)        t₂²     ]
    return Operator(ρ,(2,),true)
end

function loglikehood_gaussian(M,T,Nm)
    # Gaussian statistics for the measurement count probability
    X = t_expectation_values(T)
    return sum((Nm*p-m)^2/(2*Nm*p) for (p,m) in zip(X,M))
end

function loglikehood_binomial(M,T,Nm)
    # Binomial statistics for the measurement count probability, up to some irrelevant constant
    X = t_expectation_values(T)
    return sum(-m*log(p) for (p,m) in zip(X,M))
end

t_expectation_values(T) = (@inbounds t₁,t₂,t₃,t₄=T; t_expectation_values(t₁,t₂,t₃,t₄))

function t_expectation_values(t₁,t₂,t₃,t₄)
    # Calculate the expectation values for |+⟩,|-⟩,|+i⟩,|-i⟩,|1⟩,|0⟩ with the
    # t-parameterized density matrix ρ(t) = T†*T/trace(T†*T), where
    # T = [  t₁    0
    #      t₃+it₄  t₂]
    t₁²,t₂²,t₃²,t₄² = t₁^2,t₂^2,t₃^2,t₄^2
    N = t₁²+t₂²+t₃²+t₄²
    # ⟨+|ρ|+⟩, ⟨-|ρ|-⟩
    Xp, Xm = 1/2 + t₂*t₃/N, 1/2 - t₂*t₃/N
    # ⟨+i|ρ|+i⟩, ⟨-i|ρ|-i⟩
    Yp, Ym = 1/2 + t₂*t₄/N, 1/2 - t₂*t₄/N
    # ⟨1|ρ|1⟩, ⟨0|ρ|0⟩
    Zp, Zm = t₂²/N, (t₁²+t₃²+t₄²)/N
    return Xp,Xm,Yp,Ym,Zp,Zm
end

M = [103, 633, 621, 115, 354, 377]
