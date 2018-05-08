# Dynamics Tests
using Base.Test, Schrodinger

@testset "Time Dynamics" begin
N = 6 # cavity levels for tests

@testset "Simple Qubit" begin
ω = 0.2*2π
H = 0.5ω*σx
ψ₀ = basis(2,0)
res = sesolve(H,ψ₀,(0.0,10.0),(σz,σy),saveat=0.1)
@test real.(res.evals) ≈ [cos.(ω.*res.times) sin.(π.+ω.*res.times)]
end

@testset "Jaynes-Cummings" begin
ωc = 1.0*2π
ωa = 1.0*2π
g = 0.05*2π
κ = 0.1
a = destroy(N) ⊗ qeye(2)
sm = qeye(N) ⊗ destroy(2)
Hj = ωc*a'*a + ωa*sm'*sm + g*(a'*sm + a*sm')
O = (a'*a, sm'*sm)
ψ₀ = normalize!(basis(N,1)+0.5*basis(N,3)+1im*basis(N,4)) ⊗ basis(2,0)
res = sesolve(Hj, ψ₀, (0.0,25.0), qeye(N) ⊗ Operator([0 0;0 1]), abstol=1E-12, reltol=1E-12)
f(t,g) = 0.5*(1 - (4/9*cos(2g*t) + 1/9*cos(√3*2g*t) + 4/9*cos(√4*2g*t)))
@test real.(res.evals) ≈ f.(res.times,g)

ψ₀ = basis(N,0) ⊗ basis(2,1)
c_ops = (√(κ)*a, √(κ)*sm)
L = @inferred LindbladEvo(Hj,c_ops)
resj = lsolve(L,Operator(ψ₀),(0.0,25.0), O, Schrodinger.Tsit5(),reltol=1E-10,abstol=1E-10)
f(t,g,κ,ϕ) = 0.5*(1+cos(t*2g+ϕ))*exp(-0.5κ*t)*exp(-0.5κ*t)
@test real.(resj.evals) ≈ [f.(resj.times,g,κ,-π) f.(resj.times,g,κ,0)]
end

@testset "Three-Level Atom and Cavity" begin
c = 4.4 # coupling strength
κ = 0.1
g,e,u = basis(3,0), basis(3,1), basis(3,2)
σ_ge = qeye(N) ⊗ (g*e') # |g><e|
σ_ue = qeye(N) ⊗ (u*e') # |u><e|
a = destroy(N) ⊗ qeye(3)
n = numberop(N) ⊗ qeye(3)
ψ₀ = normalize!(basis(N,0) + 0.5*basis(N,1)) ⊗ u # Define initial state
σ_gg = Operator(basis(N,1) ⊗ g) # Define states onto which to project
σ_uu = Operator(basis(N,0) ⊗ u)
H0 = n - c * (σ_ge' * a + a' * σ_ge) # time-independent term
H1 = (σ_ue' + σ_ue);  # time-dependent term
H1_coeff(t,p) = 9*exp(-(t/5)^2)
S = @inferred SchrodingerEvo(H0,(H1,H1_coeff))
resj = lsolve(S, ψ₀, (-15.0,15.0), (n, σ_uu, σ_gg), Schrodinger.Tsit5())

S = @inferred LindbladEvo((H0,(H1,H1_coeff,[1,2,3])),√(κ)*a)
resj = lsolve(S, Operator(ψ₀), (-15.0,15.0), (n, σ_uu, σ_gg), Schrodinger.Tsit5())
end

end