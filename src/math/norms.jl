trace_norm(A::Operator) = sum(d -> √(abs(d)) , eigvals(A'*A))
