# GRAPE Amplitude Limiting via penalty
# given a set of limit control amplitudes, calculate a penalty (order ~1) when the limit is reached by the current control amplitudes. As the limit is approached, given the option for a gradually increasing penalty
#
# Questions:
# 1. given u and ulim, how do we calculate the penalty?
# 2. user function?
# 3. parameters?
# 4. do we average the penalty over all control amplitudes, or just sum them up?
#
# Discussion:
# the fidelity error used as objective function has range between 0 and 1, therefore, it makes sense for the penalty to be of order 1. (right?) We do not want the penatly to completely overtake the fidelity error, especially if we do have an option for smaller penalties near the limit.

# After much thinking, the best penalty function is:
#penalty(u,ulim=1,p=1) = p*u^10/(1000000*(ulim^10-u^10))
#penaltyprime(u,ulim=1,p=1) = 10p*abs(u)^9*ulim^10/(1000000*(u^10-ulim^10)^2)
# The parameter p is chosend such that the penalty ~ 1E-9 for u/ulim = 0.5
#penalty(u,ulim=1,p=1) = 1/(1+exp(-1000p*(u-ulim)/ulim))
#penaltyprime(u,ulim=1,p=1) = (x=penalty(u,ulim,p); 1000x*(1-x))

penalty(u,ulim=1,p=1) = (abs(u)/ulim)^100p
penaltyprime(u,ulim=1,p=1) = 100p*(abs(u)/ulim)^(100p-1)

struct AmplitudeLimiter{T<:Union{Array,Number}} <: PenaltyFunction
    limit::T
    p::Float64
end

function objective(L::AmplitudeLimiter{T},u) where {T<:Array}
    N = length(u)
    N == length(L.limit) || throw(DimensionMismatch("Control and limit amplitudes do not have the same size!"))
    res = 0.0
    for i = 1:N
        @inbounds res += penalty(u[i],L.limit[i],L.p)
    end
    return res/N
end

function objective(L::AmplitudeLimiter{T},u) where {T<:Number}
    return sum(u->penalty(u,L.limit,L.p),u)/length(u)
end

function gradient!(L::AmplitudeLimiter{T},fp,u) where {T<:Array}
    N = length(u)
    N == length(L.limit) || throw(DimensionMismatch("Control and limit amplitudes do not have the same size!"))
    for i = 1:N
        @inbounds fp[i] += penaltyprime(u[i],L.limit[i],L.p)/N
    end
    return fp
end

function gradient!(L::AmplitudeLimiter{T},fp,u) where {T<:Number}
    N = length(u)
    for i = 1:N
        @inbounds fp[i] += penaltyprime(u[i],L.limit,L.p)/N
    end
    return fp
end
