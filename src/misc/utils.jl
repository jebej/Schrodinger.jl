import Base: reverse

reverse(n::Number) = n

function gaussian(x,w)
    return exp(-4ln2*(x/w)^2)
end

function sqrtfact(n::Integer)
    if 0 <= n < 128
        return sqrtfact_table[n+1]
    elseif n >= 128 # Stirling series for values not stored in the table
        return sqrt(sqrt(2π*n)*(n/e)^n*(1.0+1/12n+1/288n^2))
    else
        throw(ArgumentError("n must be larger than or equal to 0"))
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
