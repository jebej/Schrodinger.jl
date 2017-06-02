import Base: reverse

reverse(n::Number) = n

function gaussian(x,w)
    return exp(-4ln2*(x/w)^2)
end

function sqrfact(n::Integer)
    if 0 <= n <= 20
        return sqrt(factorial(n))
    else # Stirling's formula to avoid overflow for n > 20
        return sqrt(sqrt(2π*n)*(n/e)^n)
    end
end

expim(A::AbstractMatrix) = expim!(copy(A))

function expim!(A::AbstractMatrix)
    # First decompose A into U*Λ*Uᴴ
    F = eigfact!(A)
    Λ, U = F.values, complex(F.vectors)
    # Calculate the imaginary exponential of each eigenvalue and multiply
    # each column of U to obtain B = U*exp(iΛ)
    n = length(Λ)
    B = similar(U)
    @inbounds for j = 1:n
        a = cis(Λ[j])
        @simd for i = 1:n
            B[i,j] = a * U[i,j]
        end
    end
    # Finally multiply B by Uᴴ to obtain U*exp(iΛ)*Uᴴ = exp(i*A)
    return B*U'
end

function tindex{N}(inds,sysdims::NTuple{N,Int})
    i = inds[N]
    d = 1
    @inbounds for n = reverse(1:N-1)
        d *= sysdims[n+1]
        i += d * (inds[n]-1)
    end
    return i
end
