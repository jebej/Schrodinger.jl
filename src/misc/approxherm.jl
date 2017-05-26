# The functions below check if a matrix is approximately hermitian, with a tolerance on the real and imaginary parts of HERMERR. The diagonal still must be fully real.
const HERMERR = 1E-12

function isapproxhermitian(A::AbstractMatrix)
    N,M = size(A)
    if N!==M
        return false
    end
    @inbounds for j = 1:N
        if !isreal(A[j,j]) #abs(imag(A[j,j]))>HERMERR
            return false
        end
        for i = 1:j-1
            if abs(real(A[i,j])-real(A[j,i]))>HERMERR || abs(imag(A[i,j])+imag(A[j,i]))>HERMERR
                return false
            end
        end
    end
    return true
end

function isapproxhermitian(A::SparseMatrixCSC)
    # Slighly modified from base
    m, n = size(A)
    if m != n; return false; end
    colptr = A.colptr
    rowval = A.rowval
    nzval = A.nzval
    tracker = copy(A.colptr)
    @inbounds for col = 1:A.n
        for p = tracker[col]:colptr[col+1]-1
            val = nzval[p]
            row = rowval[p]
            # Ignore stored zeros
            if val == 0
                continue
            end
            # If the matrix was symmetric we should have updated
            # the tracker to start at the diagonal or below. Here
            # we are above the diagonal so the matrix can't be symmetric.
            if row < col
                return false
            end
            # Diagonal element
            if row == col
                if !isreal(val) #abs(imag(val))>HERMERR
                    return false
                end
            else
                offset = tracker[row]
                # If the matrix is unsymmetric, there might not exist
                # a rowval[offset]
                if offset > length(rowval)
                    return false
                end
                row2 = rowval[offset]
                # row2 can be less than col if the tracker didn't
                # get updated due to stored zeros in previous elements.
                # We therefore "catch up" here while making sure that
                # the elements are actually zero.
                while row2 < col
                    if nzval[offset] != 0
                        return false
                    end
                    offset += 1
                    row2 = rowval[offset]
                    tracker[row] += 1
                end
                # Non zero A[i,j] exists but A[j,i] does not exist
                if row2 > col
                    return false
                end
                # A[i,j] and A[j,i] exists
                if row2 == col
                    if abs(real(val)-real(nzval[offset]))>HERMERR || abs(imag(val)+imag(nzval[offset]))>HERMERR
                        return false
                    end
                    tracker[row] += 1
                end
            end
        end
    end
    return true
end

function hermitianize!(A::AbstractMatrix, N=checksquare(A))
    @inbounds for j = 1:N
        A[j,j] = real(A[j,j])
         for i = 1:j-1
            A[i,j] = conj(A[j,i])
        end
    end
end
