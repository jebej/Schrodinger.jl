function gaussianpulse(t::Real,p::Vector)
    # Gaussian pulse centered on t=0, area normalized, and starting and finishing at y=0
    σ   = p[1] # gaussian standard deviation (FWHM = 2√(2ln2)σ ≈ 2.355σ)
    tg  = p[2] # total gate time
    ω10 = p[3] # qubit 0-1 transition angular freq
    ϕ   = p[4] # phase
    A   = p[5] # amplitude, (π for π-pulse, etc...)
    B = inv(√(2π)*σ*erf(tg/(√(8)*σ))-tg*gaussian(0.5tg,σ))
    Ɛπ = A*B*(gaussian(t,σ)-gaussian(0.5tg,σ))
    return Ɛπ*cos(ω10*t+ϕ)
end

function dragpulse(t::Real,p::Vector)
    # 1st order DRAG-enhanced gaussian pulse centered on t=0
    σ   = p[1] # gaussian standard deviation (FWHM = 2√(2ln2)σ ≈ 2.355σ)
    tg  = p[2] # total gate time
    ω10 = p[3] # qubit 0-1 transition angular freq
    Δ   = p[4] # ω21 - ω10
    λ   = p[5] # drive 1-2 transition amplitude ~√(2)
    ϕ   = p[6] # phase
    A   = p[7] # amplitude, (π for π-pulse, etc...)
    B = inv(√(2π)*σ*erf(tg/(√(8)*σ))-tg*gaussian(0.5tg,σ))
    Ɛπ = A*B*(gaussian(t,σ)-gaussian(0.5tg,σ))
    Ɛπ′ = A*B*gaussianprime(t,σ)
    ωt = ω10*t
    return Ɛπ*cos(ωt+ϕ) + (-λ^2*Ɛπ′/4Δ)*sin(ωt+ϕ)
    # Below is with dynamical detuning, which seems worse
    #δ = (λ^2-4)*Ɛπ^2/4Δ
    #ωt = (ω10-δ)*t
    #return Ɛπ*cos(ωt+ϕ) + (-Ɛπ′/Δ)*sin(ωt+ϕ)
end

function dragpulse5(t::Real,p::Vector)
    # 5th order DRAG-enhanced gaussian pulse centered on t=0
    σ   = p[1] # gaussian standard deviation (FWHM = 2√(2ln2)σ ≈ 2.355σ)
    tg  = p[2] # total gate time
    ω10 = p[3] # qubit 0-1 transition angular freq
    Δ   = p[4] # ω21 - ω10
    λ   = p[5] # drive 1-2 transition amplitude ~√(2)
    ϕ   = p[6] # phase
    A   = p[7] # amplitude, (π for π-pulse, etc...)
    B = inv(√(2π)*σ*erf(tg/(√(8)*σ))-tg*gaussian(0.5tg,σ))
    Ɛπ = A*B*(gaussian(t,σ)-gaussian(0.5tg,σ))
    Ɛπ′ = A*B*gaussianprime(t,σ)
    Ɛx = Ɛπ + (λ^2-4)*Ɛπ^3/8Δ^2 - (13λ^4-76λ^2+112)*Ɛπ^5/128Δ^4
    Ɛy = -λ^2*Ɛπ′/4Δ + 33(λ^2-2)*Ɛπ^2*λ^6*Ɛπ′/1536Δ^3
    ωt = ω10*t
    return Ɛx*cos(ωt+ϕ) + Ɛy*sin(ωt+ϕ)
    # Below is with dynamical detuning, which seems worse
    #Ɛy = -Ɛπ′/Δ + 33(λ^2-2)*Ɛπ^2*Ɛπ′/24Δ^3
    #δ = (λ^2-4)*Ɛπ^2/4Δ - (λ^4-7λ^2+12)*Ɛπ^4/16Δ^3
    #ωt = (ω10-δ)*t
    #return Ɛx*cos(ωt+ϕ) + Ɛy*sin(ωt+ϕ)
end


function simplegaussian(t::Real,p::Vector)
    # Gaussian pulse centered on t=0
    FWHM = p[1] # full width at half-maximum
    ω10  = p[2] # qubit 0-1 transition angular freq
    ϕ    = p[3] # phase
    A    = p[4] # amplitude
    #B = sqrtln16/(√(π)*FWHM)
    return A*gaussian_fwhm(t,FWHM)*cos(ω10*t+ϕ)
end

function simpledragpulse(t::Real,p::Vector)
    # 1st order DRAG-enhanced gaussian pulse centered on t=0
    FWHM = p[1] # full width at half-maximum
    ω10  = p[2] # qubit 0-1 transition angular freq
    Δ    = p[3] # ω21 - ω10
    λ    = p[4] # drive 1-2 transition amplitude ~√(2)
    ϕ    = p[5] # phase
    A    = p[6] # amplitude, (π for π-pulse, etc...)
    B = sqrtln16/(√(π)*FWHM)
    Ɛx = A*gaussian_fwhm(t,FWHM)*B
    Ɛy = -A*gaussianprime_fwhm(t,FWHM)*B/Δ
    δ = (λ^2-4)*Ɛx^2/4Δ
    ωt = (ω10-δ)*t
    return Ɛx*cos(ωt+ϕ) + Ɛy*sin(ωt+ϕ)
end

function simpledragpulse5(t::Real,p::Vector)
    # 5th order DRAG-enhanced gaussian pulse centered on t=0
    FWHM = p[1] # full width at half-maximum
    ω10  = p[2] # qubit 0-1 transition angular freq
    ω21  = p[3] # qubit 1-2 transition angular freq
    λ    = p[4] # drive 1-2 transition amplitude ~√(2)
    ϕ    = p[5] # phase
    A    = p[6] # amplitude, (π for π-pulse, etc...)
    Δ = ω21-ω10
    B = sqrtln16/(√(π)*FWHM)
    Ɛπ = A*gaussian_fwhm(t,FWHM)*B
    Ɛπ′ = A*gaussianprime_fwhm(t,FWHM)*B
    Ɛx = Ɛπ + (λ^2-4)*Ɛπ^3/8Δ^2 - (13λ^4-76λ^2+112)*Ɛπ^5/128Δ^4
    Ɛy = -Ɛπ′/Δ + 33(λ^2-2)*Ɛπ^2*Ɛπ′/24Δ^3
    δ = (λ^2-4)*Ɛx^2/4Δ - (λ^4-7λ^2+12)*Ɛπ^4/16Δ^3
    ωt = (ω10-δ)*t
    return Ɛx*cos(ωt+ϕ) + Ɛy*sin(ωt+ϕ)
end
