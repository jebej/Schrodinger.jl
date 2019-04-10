using SpecialFunctions

function super_basic_tomo(mx,my,mz)
    # here, the inputs are measurements in the Pauli basis
    (mx^2+my^2+mz^2)<=1 || throw(ArgumentError("invalid measurement values!"))
    ρ = [     1 + mz  mx - 1im*my
         mx + 1im*my       1 - mz]/2
    return Operator(ρ,(2,))
end


function basic_mle_tomo(M,Nm)
    # here, the inputs are e count values for |+⟩,|-⟩,|+i⟩,|-i⟩,|1⟩,|0⟩
    # meansurements, e.g. M[1] = Nm*⟨+|ρ|+⟩
    # a physical density matrix can be represented by 4 free parameters:
    # T = [  t₁    0
    #      t₃+it₄  t₂]
    # then, ρ = T†*T/trace(T†*T)
    # ρ = [t₁²+t₃²+t₄²  t₂*(t₃-it₄)    /
    #      t₂*(t₃+it₄)       t₂²   ] /   (t₁²+t₂²+t₃²+t₄²)
    # assuming gaussian noise, we can minimize the loglikehood function
    f = t -> loglikehood_gaussian(M,t[1],t[2],t[3],t[4],Nm)
    return optimize(f,rand(4))
end

function reconstruct_ρ(t₁,t₂,t₃,t₄)
    # reconstruct the density matrix from the T parameters
    t₁²,t₂²,t₃²,t₄² = t₁^2,t₂^2,t₃^2,t₄^2
    N = t₁²+t₂²+t₃²+t₄²
    ρ = 1/N*[  t₁²+t₃²+t₄²   t₂*(t₃-1im*t₄)
             t₂*(t₃+1im*t₄)        t₂²     ]
    return Operator(ρ,(2,),true)
end

function loglikehood_gaussian(M,t₁,t₂,t₃,t₄,Nm)
    # calculate the log-likelhood of measuring M, a vector of 6 expectation
    # values, given a t-parameterized density matrix
    # we assume Gaussian statistics for the measurement count probability
    X = t_expectation_values(t₁,t₂,t₃,t₄)
    return sum((Nm*px-nm)^2/(2*Nm*px) for (px,nm) in zip(X,M))
end

binomlogpdf(n::Real, p::Real, k::Real) =
    -log1p(n) - lbeta(n - k + 1, k + 1) + k * log(p) + (n - k) * log1p(-p)

function loglikehood_binomial(M,t₁,t₂,t₃,t₄,Nm)
    # calculate the log-likelhood of measuring M, a vector of 6 expectation
    # values, given a t-parameterized density matrix
    # we assume Gaussian statistics for the measurement count probability
    X = t_expectation_values(t₁,t₂,t₃,t₄)
    return sum(-binomlogpdf(Nm,px,nm) for (px,nm) in zip(X,M))
end

function t_expectation_values(t₁,t₂,t₃,t₄)
    # Calculate the expectation values for |+⟩,|-⟩,|+i⟩,|-i⟩,|1⟩,|0⟩ with the
    # t-parameterized density matrix ρ(t) = T†*T/trace(T†*T)
    # formulas calculated in Mathematica
    t₁²,t₂²,t₃²,t₄² = t₁^2,t₂^2,t₃^2,t₄^2
    N = t₁²+t₂²+t₃²+t₄²
    # |+⟩,|-⟩
    Xp, Xm = 1/2 + t₂*t₃/N, 1/2 - t₂*t₃/N
    # |+i⟩,|-i⟩
    Yp, Ym = 1/2 + t₂*t₄/N, 1/2 - t₂*t₄/N
    # |1⟩,|0⟩
    Zp, Zm = t₂²/N, (t₁²+t₃²+t₄²)/N
    return Xp,Xm,Yp,Ym,Zp,Zm
end

function
