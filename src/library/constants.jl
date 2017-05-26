export ℎ, ħ, ℏ, ℯ, Φ₀, ϕ₀, ln2,
    sigma0, σ0, sigmax, σx, sigmay, σy, sigmaz, σz, sigmaplus, σ₊, sigmaminus, σ₋,
    ⊗

# Physical constants
const ℎ = 6.626070040E-34 # Planck constant
const ħ = 1.054571800E-34 # Reduced Planck constant, in Joules*seconds, 2014 CODATA
const ℏ = ħ # Actual unicode symbol
const ℯ = 1.6021766208E-19 # Elementary charge, in Coulomb, 2014 CODATA
const Φ₀ = 2.067833831E-15 # Magnetic flux quantum (ℎ/2ℯ), in Weber, 2014 CODATA
const ϕ₀ = 3.291059757E-16 # Reduced flux quantum (ℏ/2ℯ), in Weber, 2014 CODATA
const ln2 = 0.6931471805599453 # Natural logarithm of 2

# Pauli matrices (spin-1/2 operators)
const sigma0 = operator([1.0 0.0 ; 0.0 1.0],(2,))
const sigmax = operator([0.0 1.0 ; 1.0 0.0],(2,))
const sigmay = operator([0.0 -im ; im 0.0],(2,))
const sigmaz = operator([1.0 0.0 ; 0.0 -1.0],(2,))
const σ0 = sigma0
const σx = sigmax
const σy = sigmay
const σz = sigmaz
const sigmaplus  = operator([0.0 0.0 ; 1.0 0.0],(2,)) # σ₊*[1,0] = [0,1]
const sigmaminus = operator([0.0 1.0 ; 0.0 0.0],(2,))
const σ₊ = sigmaplus
const σ₋ = sigmaminus

# Function aliases
const ⊗ = kron
const tensor = kron
