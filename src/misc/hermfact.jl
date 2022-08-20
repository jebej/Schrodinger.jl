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
