using LinearAlgebra: chkstride1, checksquare, LAPACK.chklapackerror,
    BlasInt, BLAS.@blasfunc

const liblapack = Base.liblapack_name

hermfact!(w::Vector,Z::Matrix,H::Hermitian) = syevr!(w,Z,H.uplo,H.data,0.0,0.0,0,0,-1.0)

for (syevr1, syevr2, elty, relty) in
    ((:dsyevr_,:zheevr_,:ComplexF64,:Float64),
    (:ssyevr_,:cheevr_,:ComplexF32,:Float32))
    @eval begin
        # Real symmetric
        function syevr!(w::Vector{$relty}, Z::Matrix{$relty}, uplo::AbstractChar, A::AbstractMatrix{$relty},
            vl::AbstractFloat, vu::AbstractFloat, il::Integer, iu::Integer, abstol::AbstractFloat)
            chkstride1(A)
            n = checksquare(A)
            length(w) == n || throw(ArgumentError("vector w should have length $n."))
            size(Z) == (n,n) || throw(ArgumentError("matrix Z should have size $n by $n."))
            jobz = 'V'
            range = 'A'
            lda = stride(A,2)
            m = Ref{BlasInt}()
            ldz = n
            isuppz = similar(A, BlasInt, 2*n)
            work   = Vector{$relty}(undef, 1)
            lwork  = BlasInt(-1)
            iwork  = Vector{BlasInt}(undef, 1)
            liwork = BlasInt(-1)
            info   = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1] and liwork as iwork[1]
                ccall((@blasfunc($syevr1), liblapack), Cvoid,
                (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                Ptr{$relty}, Ref{BlasInt}, Ref{$relty}, Ref{$relty},
                Ref{BlasInt}, Ref{BlasInt}, Ref{$relty}, Ptr{BlasInt},
                Ptr{$relty}, Ptr{$relty}, Ref{BlasInt}, Ptr{BlasInt},
                Ptr{$relty}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                Ptr{BlasInt}),
                jobz, range, uplo, n,
                A, max(1,lda), vl, vu,
                il, iu, abstol, m,
                w, Z, max(1,ldz), isuppz,
                work, lwork, iwork, liwork,
                info)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                    liwork = iwork[1]
                    resize!(iwork, liwork)
                end
            end
            return w, Z
        end
        # Hermitian
        function syevr!(w::Vector{$relty}, Z::Matrix{$elty}, uplo::AbstractChar, A::AbstractMatrix{$elty},
            vl::AbstractFloat, vu::AbstractFloat, il::Integer, iu::Integer, abstol::AbstractFloat)
            chkstride1(A)
            n = checksquare(A)
            length(w) == n || throw(ArgumentError("vector w should have length $n."))
            size(Z) == (n,n) || throw(ArgumentError("matrix Z should have size $n by $n."))
            jobz = 'V'
            range = 'A'
            lda = max(1,stride(A,2))
            m = Ref{BlasInt}()
            ldz = n
            isuppz = similar(A, BlasInt, 2*n)
            work   = Vector{$elty}(undef, 1)
            lwork  = BlasInt(-1)
            rwork  = Vector{$relty}(undef, 1)
            lrwork = BlasInt(-1)
            iwork  = Vector{BlasInt}(undef, 1)
            liwork = BlasInt(-1)
            info   = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1], lrwork as rwork[1] and liwork as iwork[1]
                ccall((@blasfunc($syevr2), liblapack), Cvoid,
                (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                Ptr{$elty}, Ref{BlasInt}, Ref{$elty}, Ref{$elty},
                Ref{BlasInt}, Ref{BlasInt}, Ref{$elty}, Ptr{BlasInt},
                Ptr{$relty}, Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt},
                Ptr{$elty}, Ref{BlasInt}, Ptr{$relty}, Ref{BlasInt},
                Ptr{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}),
                jobz, range, uplo, n,
                A, lda, vl, vu,
                il, iu, abstol, m,
                w, Z, ldz, isuppz,
                work, lwork, rwork, lrwork,
                iwork, liwork, info)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                    lrwork = BlasInt(rwork[1])
                    resize!(rwork, lrwork)
                    liwork = iwork[1]
                    resize!(iwork, liwork)
                end
            end
            return w, Z
        end
    end
end
