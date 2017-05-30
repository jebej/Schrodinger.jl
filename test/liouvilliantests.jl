using Schrodinger, BenchmarkTools

N = 10   # Set where to truncate Fock state for cavity
g = 4.4; # coupling strength
gs = basis(3, 0)
es = basis(3, 1)
us  = basis(3, 2)
σ_ge = tensor(qeye(N), gs*es')  # |g><e|
σ_ue = tensor(qeye(N), us*es')  # |u><e|
a = tensor(destroy(N), qeye(3))
ada = tensor(numberop(N), qeye(3))
psi0 = tensor(normalize!(basis(N, 0)+0.5*basis(N, 1)), us) # Define initial state
σ_gg = Operator(tensor(basis(N, 1), gs)) # Define states onto which to project
σ_uu = Operator(tensor(basis(N, 0), us))
H0 = ada - g * (σ_ge' * a + a' * σ_ge)  # time-independent term
H1 = (σ_ue' + σ_ue)  # time-dependent term

H1_coeff(t,params)= 9 * exp(-(t / 5.)^ 2)

L = Schrodinger.SchrodingerEvo([H0,(H1,H1_coeff)])


t  = 0.1
dψ = similar(full(psi0))
ψ  = copy(full(psi0))

L(t,ψ,dψ)
