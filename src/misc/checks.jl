# The functions below check if a matrix is approximately Hermitian
const ERR = 1E-13

function isapproxhermitian(A::AbstractMatrix)
    m, n = size(A)
    m == n || return false
    @inbounds for j = 1:m
        for i = 1:j-1
            Aij, Aji = A[i,j], A[j,i]
            abs(real(Aij)-real(Aji)) > ERR || abs(imag(Aij)+imag(Aji)) > ERR && return false
        end
        abs(imag(A[j,j])) > ERR && return false
    end
    return true
end

function isapproxhermitian(A::SparseMatrixCSC)
    # Slighly modified from base
    m, n = size(A)
    m == n || return false
    colptr = A.colptr
    rowval = A.rowval
    nzval = A.nzval
    tracker = copy(A.colptr)
    @inbounds for col = 1:m
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
                abs(imag(val)) > ERR && return false
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
                    abs(real(val)-real(nzval[offset])) > ERR || abs(imag(val)+imag(nzval[offset])) > ERR && return false
                    tracker[row] += 1
                end
            end
        end
    end
    return true
end

function hermitianize!(A::AbstractMatrix)
    m = checksquare(A)
    @inbounds for j = 1:m
        A[j,j] = real(A[j,j])
         for i = 1:j-1
            A[i,j] = conj(A[j,i])
        end
    end
end
