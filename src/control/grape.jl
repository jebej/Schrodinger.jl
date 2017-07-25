function hsinner!(U1::Operator,U2::Operator,C)
    Ac_mul_B!(C,data(U1),data(U2))
    return trace(C)/prod(dims(U1))
end

function calc_prop_steps(Hd::Operator,Hc::Vector{Operator},u::Matrix{Real},t::Real,n::Integer)
    m = length(Hc)
    m,n == size(u) || throw(ArgumentError("control amplitude matrix not consistent with number of timesteps or number of control hamiltonians"))
    N = prod(dims(Hd))
    δt = t/n
    Hj = Matrix{Complex128}(N,N)
    UU = map(1:n) do j
        copy!(Hj,full(Hd))
        for k = 1:m
            Hj .+= u[k,j].*full(Hc[k])
        end
        return expim!(scale!(Hj,-δt))
    end
end
