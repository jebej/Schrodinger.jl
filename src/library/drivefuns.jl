function gaussianpulse(t::Real,p::Vector)
    # Gaussian pulse centered on t=0, area normalized, and starting and finishing at y=0
    σ   = p[1] # gaussian standard deviation (FWHM = 2√(2ln2)σ ≈ 2.355σ)
    tg  = p[2] # total gate time
    ω10 = p[3] # qubit 0-1 transition angular freq
    ϕ   = p[4] # phase
    A   = p[5] # amplitude, (π for π-pulse, etc...)
    B = inv(√(2π)*σ*erf(tg/(√(8)*σ))-tg*gaussian(tg/2,σ))
    Ɛπ = A*B*(gaussian(t,σ)-gaussian(tg/2,σ))
    return Ɛπ*cos(ω10*t+ϕ)
end

function cosinepulse(t::Real,p::Vector)
    # Cosine pulse centered on t=0
    σ   = p[1] # not used
    tg  = p[2] # total gate time
    ω10 = p[3] # qubit 0-1 transition angular freq
    ϕ   = p[4] # phase
    A   = p[5] # amplitude, (π for π-pulse, etc...)
    Ɛπ = A/tg*(cos(2π*t/tg)+1)
    return Ɛπ*cos(ω10*t+ϕ)
end

function dragpulse(t::Real,p::Vector)
    # 2nd order DRAG-enhanced Gaussian pulse (Y-only correction) centered on t=0
    # 10.1103/PhysRevA.83.012308
    σ   = p[1] # gaussian standard deviation (FWHM = 2√(2ln2)σ ≈ 2.355σ)
    tg  = p[2] # total gate time
    ω10 = p[3] # qubit 0-1 transition angular freq
    Δ   = p[4] # ω21 - ω10
    λ   = p[5] # drive 1-2 transition amplitude ~√(2)
    ϕ   = p[6] # phase
    A   = p[7] # amplitude, (π for π-pulse, etc...)
    B = inv(√(2π)*σ*erf(tg/(√(8)*σ))-tg*gaussian(tg/2,σ))
    Ɛπ = A*B*(gaussian(t,σ)-gaussian(tg/2,σ))
    Ɛπ′ = A*B*gaussianprime(t,σ)
    ωt = ω10*t
    return (Ɛπ - λ^2*(λ^2-4)*Ɛπ^3/32Δ^2)*cos(ωt+ϕ) + (-λ^2*Ɛπ′/4Δ)*sin(ωt+ϕ)
end

function dragpulse_det(t::Real,p::Vector)
    # 2nd order DRAG-enhanced Gaussian pulse (Y-only) centered on t=0
    # with detuning
    # PhysRevLett.116.020501
    σ   = p[1] # gaussian standard deviation (FWHM = 2√(2ln2)σ ≈ 2.355σ)
    tg  = p[2] # total gate time
    ω10 = p[3] # qubit 0-1 transition angular freq
    Δ   = p[4] # ω21 - ω10
    λ   = p[5] # drive 1-2 transition amplitude ~√(2)
    δω  = p[6]
    ϕ   = p[7] # phase
    A   = p[8] # amplitude, (π for π-pulse, etc...)
    Δ = Δ - δω
    B = inv(√(2π)*σ*erf(tg/(√(8)*σ))-tg*gaussian(tg/2,σ))
    Ɛπ = A*B*(gaussian(t,σ)-gaussian(tg/2,σ))
    Ɛπ′ = A*B*gaussianprime(t,σ)
    C = (Ɛπ - λ^2*(λ^2-4)*Ɛπ^3/32Δ^2)
    S = (-λ^2*Ɛπ′/4Δ)
    I = C*cos(δω*t) - S*sin(δω*t)
    Q = S*cos(δω*t) + C*sin(δω*t)
    ωt = ω10*t
    return I*cos(ωt+ϕ) + Q*sin(ωt+ϕ)
end

function dragpulse_cos(t::Real,p::Vector)
    # 2nd order DRAG-enhanced cosine pulse (Y-only) centered on t=0
    # 10.1103/PhysRevA.83.012308
    σ   = p[1] # gaussian standard deviation (FWHM = 2√(2ln2)σ ≈ 2.355σ)
    tg  = p[2] # total gate time
    ω10 = p[3] # qubit 0-1 transition angular freq
    Δ   = p[4] # ω21 - ω10
    λ   = p[5] # drive 1-2 transition amplitude ~√(2)
    ϕ   = p[6] # phase
    A   = p[7] # amplitude, (π for π-pulse, etc...)
    B = 1
    Ɛπ = A*(cos(2π*t/tg)+1)/tg
    Ɛπ′ = -A*sin(2π*t/tg)*2π/tg^2
    ωt = ω10*t
    return (Ɛπ - λ^2*(λ^2-4)*Ɛπ^3/32Δ^2)*cos(ωt+ϕ) + (-λ^2*Ɛπ′/4Δ)*sin(ωt+ϕ)
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
    B = inv(√(2π)*σ*erf(tg/(√(8)*σ))-tg*gaussian(tg/2,σ))
    Ɛπ = A*B*(gaussian(t,σ)-gaussian(tg/2,σ))
    Ɛπ′ = A*B*gaussianprime(t,σ)
    Ɛx = Ɛπ + (λ^2-4)*Ɛπ^3/8Δ^2 - (13λ^4-76λ^2+112)*Ɛπ^5/128Δ^4
    Ɛy = -Ɛπ′/Δ + 33(λ^2-2)*Ɛπ^2*Ɛπ′/24Δ^3
    δ = (λ^2-4)*Ɛπ^2/4Δ - (λ^4-7λ^2+12)*Ɛπ^4/16Δ^3
    ωt = (ω10-δ)*t
    return Ɛx*cos(ωt+ϕ) + Ɛy*sin(ωt+ϕ)
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

function stepfun(t::Real,p::Vector,u::Vector)
    tg = p[1] # total gate time
    to = p[2] # time offset (to make the time start at 0, e.g. if tspan=(-2,3), to=2)
    return u[round(Int,1+(length(u)-1)*(t+to)/tg)]
end
