using FastLapackInterface

const HERMEIGENWSREAL = Ref{HermitianEigenWs{Float64,Matrix{Float64},Float64}}()
const HERMEIGENWSCOMP = Ref{HermitianEigenWs{ComplexF64,Matrix{ComplexF64},Float64}}()

function _get_syevr_ws(A::AbstractMatrix{T}) where T<:Real
    isdefined(HERMEIGENWSREAL,1) || setindex!(HERMEIGENWSREAL, HermitianEigenWs(A, vecs=true))
    return HERMEIGENWSREAL[]
end
function _get_syevr_ws(A::AbstractMatrix{T}) where T<:Complex
    isdefined(HERMEIGENWSCOMP,1) || setindex!(HERMEIGENWSCOMP, HermitianEigenWs(A, vecs=true))
    return HERMEIGENWSCOMP[]
end

hermfact!(H::Hermitian) = LAPACK.syevr!(_get_syevr_ws(H.data),'V','A',H.uplo,H.data,0.0,0.0,0,0,-1.0,resize=true)

# Change of basis to convert from a Hermitian Tridiagonal matrix to a real symmetric Tridiagonal matrix
# ref: https://groups.google.com/g/comp.lang.fortran/c/UF1D1D27xJQ/m/j_nHpB3lsdEJ

hermitian_tridiag_to_realsym_tridiag(H::Tridiagonal) = hermitian_tridiag_to_realsym_tridiag(H.d, H.dl)

function hermitian_tridiag_to_realsym_tridiag(dv::AbstractVector, dl::AbstractVector)
    # dv: diagonal values, dl: subdiagonal values
    # no check is done to verify H = Tridiagonal(dl,dv,conj.dl)) is indeed Hermitian
    # dv should be real, the superdiagonal values should be du = conj(dl)
    # then, T = D'*H*D, where D is a diagonal (unitary) matrix and T
    # is a real symmetrix tridiagonal matrix
    N = length(dv)
    D = similar(dl, N)
    D[1] = 1
    @inbounds for k = 2:N
        if dl[k-1] == 0
            D[k] = D[k-1]
        else
            D[k] = D[k-1] * dl[k-1] / abs(dl[k-1])
        end
    end
    T = SymTridiagonal(real(dv), abs.(dl))
    return T, Diagonal(D)
end
