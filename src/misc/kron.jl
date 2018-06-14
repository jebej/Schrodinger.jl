function I_kron_A{T}(d::Int, A::AbstractMatrix{T})
    n, m = size(A)
    B = zeros(T,d*n,d*m)
    @inbounds for k = 1:d
        i = (k-1)*n + 1
        j = (k-1)*m + 1
        B[i:i+n-1, j:j+m-1] .= A
    end
    return B
end

function I_kron_A_mul_B!(C,A,B)
    n = checksquare(A)
    n^2 == checksquare(B) == checksquare(C) || throw(DimensionMismatch("invalid dimensions"))
    for j = 1:n
        jj = (j-1)*n+1 : (j-1)*n+n
        for i = 1:n
            ii = (i-1)*n+1 : (i-1)*n+n
            Cij = view(C,ii,jj)
            Bij = view(B,ii,jj)
            A_mul_B!(Cij,A,Bij)
        end
    end
    return C
end

# Not too useful:
function At_kron_I_mul_B!(C,A,B)
    n = checksquare(A)
    n^2 == checksquare(B) == checksquare(C) || throw(DimensionMismatch("invalid dimensions"))
    fill!(C,0)
    kk = 1:n
    for j = 1:n
        jj = (j-1)*n+1 : (j-1)*n+n
        for i = 1:n
            ii = (i-1)*n+1 : (i-1)*n+n
            Cij = view(C,ii,jj)
            for k = 1:n
                Cij .+= A[k,i] .* view(B,kk+(k-1)*n,jj)
            end
        end
    end
    return C
end
