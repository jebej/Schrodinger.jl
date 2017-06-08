function expect{T,D}(σ::Operator{NonHerm,T,D},ψ::Ket)
    dimsmatch(σ,ψ)
    return dot(ψ.data,σ.data*ψ.data)
end
function expect{T,D}(σ::Operator{Herm,T,D},ψ::Ket)
    dimsmatch(σ,ψ)
    return real(dot(ψ.data,σ.data*ψ.data))
end
expect(ψ::Ket,σ::Operator) = expect(σ,ψ)

function expect{T1,T2,D}(σ::Operator{NonHerm,T1,D},ρ::Operator{Herm,T2,D})
    dimsmatch(σ,ρ)
    return trace(σ.data*ρ.data)
end

function expect{T1,T2,D}(σ::Operator{Herm,T1,D},ρ::Operator{NonHerm,T2,D})
    dimsmatch(σ,ρ)
    return trace(σ.data*ρ.data)
end

function expect{T1,T2,D}(σ::Operator{Herm,T1,D},ρ::Operator{Herm,T2,D})
    dimsmatch(σ,ρ)
    return real(trace(σ.data*ρ.data))
end
