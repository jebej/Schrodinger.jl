# The functions below check if a matrix is approximately Hermitian
using Compat.LinearAlgebra: RealHermSymComplexHerm

const ERR = 1E-13

isapproxhermitian(A::RealHermSymComplexHerm) = true

function isapproxhermitian(A::AbstractMatrix)
    m, n = size(A)
    m == n || return false
    @inbounds for j = 1:m
        for i = 1:j-1
            Aij, Aji = A[i,j], A[j,i]
            if abs(real(Aij)-real(Aji)) > ERR || abs(imag(Aij)+imag(Aji)) > ERR
                return false
            end
        end
        abs(imag(A[j,j])) > ERR && return false
    end
    return true
end

function isapproxhermitian(A::SparseMatrixCSC)
    # Modified from base
    m, n = size(A)
    m == n || return false
    colptr = A.colptr
    rowval = A.rowval
    nzval = A.nzval
    @inbounds for col = 1:m
        for p = colptr[col]:colptr[col+1]-1
            val = nzval[p]
            row = rowval[p]
            if val == 0 # ignore stored zeros
                continue
            end
            if row == col # diagonal element
                abs(imag(val)) > ERR && return false
            else # off-diagonal element
                val2 = A[col,row]
                if abs(real(val)-real(val2)) > ERR || abs(imag(val)+imag(val2)) > ERR
                    return false
                end
            end
        end
    end
    return true
end

function hermitianize!(A::AbstractMatrix)
    m = checksquare(A)
    @inbounds for j = 1:m
         for i = 1:j-1
            A[i,j] = conj(A[j,i])
        end
        A[j,j] = real(A[j,j])
    end
    return A
end

function isunitary(A::AbstractMatrix{T}) where {T}
    # check efficiently if A'*A = I
    m, n = size(A)
    m == n || return false
    @inbounds for i = 1:m, j = 1:m
        a = zero(T)
        for k = 1:m
            a += conj(A[k,i]) * A[k,j]
        end
        if i == j
            a == 1 || return false
        else
            a == 0 || return false
        end
    end
    return true
end

function isunitary(A::SparseMatrixCSC{T,S}) where {T,S}
    # check efficiently if A'*A = I
    m, n = size(A)
    m == n || return false
    colptr = A.colptr
    rowval = A.rowval
    nzval = A.nzval
    @inbounds for i = 1:m, j = 1:m
        a = zero(T)
        for p1 = colptr[i]:colptr[i+1]-1, p2 = colptr[j]:colptr[j+1]-1
            if rowval[p1] < rowval[p2]
                continue
            elseif rowval[p1]==rowval[p2]
                a += conj(nzval[p1]) * nzval[p2]
                continue
            end
        end
        if i == j
            a == 1 || return false
        else
            a == 0 || return false
        end
    end
    return true
end

function isapproxunitary(A::AbstractMatrix{T}) where {T}
    # check efficiently if A'*A ≈ I
    m, n = size(A)
    m == n || return false
    @inbounds for i = 1:m, j = 1:m
        a = zero(T)
        for k = 1:m
            a += conj(A[k,i]) * A[k,j]
        end
        if i == j
            if abs(1-real(a)) > ERR || abs(imag(a)) > ERR
                return false
            end
        else
            if abs(real(a)) > ERR || abs(imag(a)) > ERR
                return false
            end
        end
    end
    return true
end

function isapproxunitary(A::SparseMatrixCSC{T,S}) where {T,S}
    # check efficiently if A'*A ≈ I
    m, n = size(A)
    m == n || return false
    colptr = A.colptr
    rowval = A.rowval
    nzval = A.nzval
    @inbounds for i = 1:m, j = 1:m
        a = zero(T)
        for p1 = colptr[i]:colptr[i+1]-1, p2 = colptr[j]:colptr[j+1]-1
            if rowval[p1] < rowval[p2]
                continue
            elseif rowval[p1]==rowval[p2]
                a += conj(nzval[p1]) * nzval[p2]
                continue
            end
        end
        if i == j
            if abs(1-real(a)) > ERR || abs(imag(a)) > ERR
                return false
            end
        else
            if abs(real(a)) > ERR || abs(imag(a)) > ERR
                return false
            end
        end
    end
    return true
end
