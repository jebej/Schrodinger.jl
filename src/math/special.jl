function expect(σ::Operator,ψ::Ket)
    dimsmatch(σ,ψ)
    val = dot(ψ.data,σ.data*ψ.data)
    #return ishermitian(σ) ? real(val) : val
end
expect(ψ::Ket,σ::Operator) = expect(σ,ψ)

function expect(σ::Operator,ρ::Density)
    dimsmatch(σ,ρ)
    val = trace(σ.data*ρ.data)
    #return ishermitian(σ) ? real(val) : val
end
expect(ρ::Density,σ::Operator) = expect(σ,ρ)
