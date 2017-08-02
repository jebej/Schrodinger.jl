using Base.LinAlg: chkstride1, checksquare, LAPACK.chklapackerror
const liblapack = Base.liblapack_name

hermfact!(w::Vector,Z::Matrix,H::Hermitian) = syevr!(w,Z,H.uplo,H.data,0.0,0.0,0,0,-1.0)

function syevr!(w::Vector{Float64}, Z::Matrix{Complex128}, uplo::Char, A::StridedMatrix{Complex128},
                vl::AbstractFloat, vu::AbstractFloat, il::Integer, iu::Integer, abstol::AbstractFloat)
    chkstride1(A)
    n = checksquare(A)
    length(w) == n || throw(ArgumentError("vector w should have length $n."))
    size(Z) == (n,n) || throw(ArgumentError("matrix Z should have size $n by $n."))
    jobz = 'V'
    range = 'A'
    lda = max(1,stride(A,2))
    m = Ref{Int}()
    ldz = n
    isuppz = Vector{Int}(2n)
    work   = Vector{Complex128}(1)
    lwork  = Int(-1)
    rwork  = Vector{Float64}(1)
    lrwork = Int(-1)
    iwork  = Vector{Int}(1)
    liwork = Int(-1)
    info   = Ref{Int}()
    for i = 1:2
        ccall((:zheevr_64_, liblapack), Void,
              (Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, Ptr{Int},
               Ptr{Complex128}, Ptr{Int}, Ptr{Complex128}, Ptr{Complex128},
               Ptr{Int}, Ptr{Int}, Ptr{Complex128}, Ptr{Int},
               Ptr{Float64}, Ptr{Complex128}, Ptr{Int}, Ptr{Int},
               Ptr{Complex128}, Ptr{Int}, Ptr{Float64}, Ptr{Int},
               Ptr{Int}, Ptr{Int}, Ptr{Int}),
              &jobz, &range, &uplo, &n,
              A, &lda, &vl, &vu,
              &il, &iu, &abstol, m,
              w, Z, &ldz, isuppz,
              work, &lwork, rwork, &lrwork,
              iwork, &liwork, info)
        chklapackerror(info[])
        if i == 1
            lwork = Int(real(work[1]))
            work = Vector{Complex128}(lwork)
            lrwork = Int(rwork[1])
            rwork = Vector{Float64}(lrwork)
            liwork = iwork[1]
            iwork = Vector{Int}(liwork)
        end
    end
    return w, Z
end

function syevr!(w::Vector{Float64}, Z::Matrix{Float64}, uplo::Char, A::StridedMatrix{Float64},
                vl::AbstractFloat, vu::AbstractFloat, il::Integer, iu::Integer, abstol::AbstractFloat)
    chkstride1(A)
    n = checksquare(A)
    length(w) == n || throw(ArgumentError("vector w should have length $n."))
    size(Z) == (n,n) || throw(ArgumentError("matrix Z should have size $n by $n."))
    jobz = 'V'
    range = 'A'
    lda = max(1,stride(A,2))
    m = Ref{Int}()
    ldz = n
    isuppz = Vector{Int}(2n)
    work   = Vector{Float64}(1)
    lwork  = Int(-1)
    iwork  = Vector{Int}(1)
    liwork = Int(-1)
    info   = Ref{Int}()
    for i = 1:2
        ccall((:dsyevr_64_, liblapack), Void,
            (Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, Ptr{Int},
                Ptr{Float64}, Ptr{Int}, Ptr{Float64}, Ptr{Float64},
                Ptr{Int}, Ptr{Int}, Ptr{Float64}, Ptr{Int},
                Ptr{Float64}, Ptr{Float64}, Ptr{Int}, Ptr{Int},
                Ptr{Float64}, Ptr{Int}, Ptr{Int}, Ptr{Int},
                Ptr{Int}),
            &jobz, &range, &uplo, &n,
            A, &lda, &vl, &vu,
            &il, &iu, &abstol, m,
            w, Z, &ldz, isuppz,
            work, &lwork, iwork, &liwork,
            info)
        chklapackerror(info[])
        if i == 1
            lwork = Int(real(work[1]))
            work = Vector{Float64}(lwork)
            liwork = iwork[1]
            iwork = Vector{Int}(liwork)
        end
    end
    return w, Z
end
