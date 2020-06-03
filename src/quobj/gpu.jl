gpu(A::Operator, GPUType) = Operator(GPUType(data(A)),dims(A))

gpu(A::Ket, GPUType) = Ket(GPUType(data(A)),dims(A))

gpu(L::Liouvillian, GPUType) = Liouvillian(L.dims, GPUType(SparseMatrixCSC{ComplexF32}(L.L₀)), GPUType.(SparseMatrixCSC{ComplexF32}.(L.Lₙ)), L.fₙ, L.pₙ)
