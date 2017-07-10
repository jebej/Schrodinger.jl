scale!(A::LinAlg.HermOrSym,b::Number) = (scale!(A.data,b);A)
