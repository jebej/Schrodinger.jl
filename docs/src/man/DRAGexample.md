```@meta
DocTestSetup  = quote
    using Schrodinger
end
```
```@setup plots
using Schrodinger, PyPlot

function rotgaussianpulse(t::Real,p::Vector)
    # normalized pulse centered on t=0, begins and ends at 0
    tg  = p[1] # gate time
    σ   = p[2] # standard deviation
    A   = p[3] # amplitude
    B   = inv(√(2π)*σ*erf(tg/(√(8)*σ))-tg*gaussian(0.5tg,σ))
    Ɛˣ = A*B*(gaussian(t,σ)-gaussian(0.5tg,σ))
    return Ɛˣ
end

function DRAGx(t::Real,p::Vector)
    # Ɛˣ for fifth order DRAG
    tg  = p[1] # gate time
    σ   = p[2] # standard deviation
    A   = p[3] # amplitude
    Δ   = p[4] # anharmonicity
    λ   = p[5] # relative strength of transitions
    Ɛˣ = rotgaussianpulse(t,[tg,σ,A]) + (λ^2-4)*rotgaussianpulse(t,[tg,σ,A])^3/(8*Δ^2) - (13λ^4-76λ^2+112)*rotgaussianpulse(t,[tg,σ,A])^5/(128Δ^4)
    return Ɛˣ
end

function DRAGy(t::Real,p::Vector)
    # Ɛʸ for fifth order DRAG
    tg  = p[1] # gate time
    σ   = p[2] # standard deviation
    A   = p[3] # amplitude
    Δ   = p[4] # anharmonicity
    λ   = p[5] # relative strength of transitions
    B   = inv(√(2π)*σ*erf(tg/(√(8)*σ))-tg*gaussian(0.5tg,σ))
    Ɛˣ′ = A*B*Schrodinger.gaussianprime(t,σ)
    Ɛʸ = -Ɛˣ′/Δ + 33*(λ^2-2)*rotgaussianpulse(t,[tg,σ,A])^2*Ɛˣ′/(24*Δ^3)
    return Ɛʸ
end

function detuning(t::Real,p::Vector)
    # dynamical detuning for fifth order DRAG
    tg  = p[1] # gate time
    σ   = p[2] # standard deviation
    A   = p[3] # amplitude
    Δ   = p[4] # anharmonicity
    λ   = p[5] # relative strength of transitions
    δ₁ = (λ^2-4)*rotgaussianpulse(t,[tg,σ,A])^2/(4*Δ) - (λ^4-7λ^2+12)*rotgaussianpulse(t,[tg,σ,A])^4/(16Δ^3)
    return δ₁
end

H = σx/2 # same as create(2)/2 + destroy(2)/2; using natural units where ħ=1

g = basis(2,0) # begin in ground state
tg = 6e-9 # gate time
σ = 0.5tg # standard deviation of gaussian
tspan = (-tg/2,tg/2) # pulse centred at t=0

res1 = sesolve([qzero(2),(H,rotgaussianpulse,[tg,σ,π])],g,tspan,saveat=linspace(-tg/2,tg/2,100))
!isdir("img") && mkdir("img")
plot(res1.times*1e9,levelprobs(res1.states)); xlabel("Time (ns)"); ylabel("Level Probabilities"); legend(["Ground State", "Excited State"]); grid()
savefig(joinpath("img","qubitNOT.svg"))
clf()

Δ = 2π*(-400e6) # anharmonicity
λ = √2 # relative transition strength
Π₂ = basis(3,2) * basis(3,2)' # projector for the 2nd level
Hc = Δ*Π₂ # constant Hamiltonian
Hd = create(3)/2 + destroy(3)/2 # affected by Ɛˣ
g = basis(3,0)

res2 = sesolve([Hc,(Hd,rotgaussianpulse,[tg,σ,π])],g,tspan)
plot(res2.times*1e9,levelprobs(res2.states)); xlabel("Time (ns)"); ylabel("Level Probabilities"); legend(["Ground State", "1st Excited State", "2nd Excited State"]); grid()
savefig(joinpath("img","3levelNOT.svg"))
clf()

Π₁ = basis(3,1) * basis(3,1)' # projector for 1st level
Hdet = Π₁ # affected by dynamical detuning
Hdx = create(3)/2 + destroy(3)/2 # affected by Ɛˣ
Hdy = im*create(3)/2 - im*destroy(3)/2 # affected by Ɛʸ

res3 = sesolve([Hc,(Hdet,detuning,[tg,σ,π,Δ,λ]),(Hdx,DRAGx,[tg,σ,π,Δ,λ]),(Hdy,DRAGy,[tg,σ,π,Δ,λ])],g,tspan)
plot(res3.times*1e9,levelprobs(res3.states)); xlabel("Time (ns)"); ylabel("Level Probabilities"); legend(["Ground State", "1st Excited State", "2nd Excited State"]); grid()
savefig(joinpath("img","3levelDRAG.svg"))
clf()

tgs = (3:9)*1e-9 # gate times
Fg_res = Matrix{Float64}(length(tgs),2) # initialize matrix for number of solves
axialStates  = [normalize!(Ket([1,1,0])),   # +X
                normalize!(Ket([1,-1,0])),  # -X
                normalize!(Ket([1,im,0])),  # +Y
                normalize!(Ket([1,-im,0])), # -Y
                normalize!(Ket([1,0,0])),   # +Z
                normalize!(Ket([0,1,0]))]   # -Z
axialOperators = [] # need density operators too
for i = 1:6
    push!(axialOperators,axialStates[i]*axialStates[i]')
end
Uideal = complex(qzero(3)); Uideal[1:2,1:2] = data(σx)
for (i,tg) in enumerate(tgs)
    sum1 = 0 # sum of Gaussian gate fidelities
    sum2 = 0 # sum of DRAG gate fidelities
    for j = 1:6
        # Gaussian
        res4_1 = sesolve([Hc,(Hd,rotgaussianpulse,[tg,σ,π])],axialStates[j],(-tg/2,tg/2))
        sum1 += trace(Uideal*axialOperators[j]*Uideal'*(res4_1.states[end]*res4_1.states[end]'))
        # DRAG
        res4_2 = sesolve([Hc,(Hdet,detuning,[tg,σ,π,Δ,λ]),(Hdx,DRAGx,[tg,σ,π,Δ,λ]),(Hdy,DRAGy,[tg,σ,π,Δ,λ])],axialStates[j],(-tg/2,tg/2))
        sum2 += trace(Uideal*axialOperators[j]*Uideal'*(res4_2.states[end]*res4_2.states[end]'))
    end
    Fg_res[i,:] = [sum1/6 sum2/6] # take average
end
plot(tgs*1e9,1.-Fg_res); ylim([10e-8,1]); title("Average gate fidelity averaging over all input states"); yscale("log"); xlabel("Gate Time (ns)"); ylabel("Gate Error 1-Fg"); legend(["Gaussian","DRAG 5th Order"]); grid()
savefig(joinpath("img","fidelities.svg"))
```

# DRAG

This example shows how to implement a NOT gate on a qubit using a Gaussian pulse. We then extend this method to a three-level slightly anharmonic energy spectrum with nearest level coupling. As we will see, the gate error increases due to leakage into the third level. To remedy this, we implement Derivative Removal by Adiabatic Gate (DRAG) which offers better gate fidelity than an ordinary Gaussian pulse \[[1]].

## A simple NOT Gate

First, we will apply a simple NOT gate to a qubit in the ground state. The hamiltonian for our qubit in the lab frame can be written as $$ħ(ω|1⟩⟨1|+ℇ(t)σ<sub>x</sub>)$$ where $$ħω$$ is the transition energy, $$σ<sub>x</sub>$$ is the Pauli-X operator, and $$ℇ(t)=ℇ<sup>x</sup>(t)cos(ω<sub>d</sub>t)$$ represents our control of the system using a drive frequency $$ω<sub>d</sub>$$. Any control $$ℇ<sup>x</sup>(t)$$ such that the integral of $$ℇ<sup>x</sup>$$ over the total gate time equals $$π$$ will result in a complete inversion. It is common to use a Gaussian shaped π-pulse to implement a NOT gate. For our purposes, we will find it more convenient to work in the rotating frame with respect to our drive frequency $$ω<sub>d</sub>$$. When this frequency is resonant with the qubit frequency $$ω$$, the Hamiltonian is given by $$ħℇ<sup>x</sup>(t)σ<sub>x</sub>/2$$. If we are going to solve the time dynamics for this system, we also have to define our pulse in the rotating frame.
```jldoctest example1
function rotgaussianpulse(t::Real,p::Vector)
    # normalized pulse centered on t=0, begins and ends at 0
    tg  = p[1] # gate time
    σ   = p[2] # standard deviation
    A   = p[3] # amplitude
    B   = inv(√(2π)*σ*erf(tg/(√(8)*σ))-tg*gaussian(0.5tg,σ)) # normalize
    Ɛˣ = A*B*(gaussian(t,σ)-gaussian(0.5tg,σ))
    return Ɛˣ
end

H = σx/2 # same as create(2)/2 + destroy(2)/2; using natural units where ħ=1
# output
2×2 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with space dimensions 2:
 0.0  0.5
 0.5  0.0
```

We are now ready to solve for the time dynamics and plot our results. We'll use parameters that are typical for Transmon qubits.

```jldoctest example1
tg = 6e-9 # gate time 6ns
σ = 0.5tg # standard deviation of gaussian pulse
g = basis(2,0) # begin in ground state
tspan = (-tg/2,tg/2) # pulse centred at t=0

res1 = sesolve([qzero(2),(H,rotgaussianpulse,[tg,σ,π])],g,tspan,saveat=linspace(-tg/2,tg/2,100))
plot(res1.times*1e9,levelprobs(res1.states)); xlabel("Time (ns)"); ylabel("Level Probabilities"); legend(["Ground State", "Excited State"]); grid();
```
![qubit NOT gate](img/qubitNOT.svg)

As is expected, the system moves from the ground state to the excited state.

## Three-level system

In reality, one must worry about leakage into other states of the system, especially when dealing with short gate times. Let's see how our Gaussian pulse performs on a three-level anharmonic system. We will again work in the rotating frame at resonance with the qubit frequency. This system will have an anharmonicity $$Δ$$, which is the detuning of the 2nd excited state with respect to the drive frequency, and an additional parameter $$λ$$ describing the relative strength of the 1-2 transition compared to the 0-1 transition (see \[[1]] for more details). Let's create the Hamiltonian and perform the same time evolution.

```jldoctest example1
Δ = 2π*(-400e6) # anharmonicity
λ = √2 # relative transition strength
Π₂ = basis(3,2) * basis(3,2)' # projector for the 2nd level
Hc = Δ*Π₂ # constant Hamiltonian
Hd = create(3)/2 + destroy(3)/2 # affected by Ɛˣ
g = basis(3,0)

res2 = sesolve([Hc,(Hd,rotgaussianpulse,[tg,σ,π])],g,tspan)
figure()
plot(res2.times*1e9,levelprobs(res2.states)); xlabel("Time (ns)"); ylabel("Level Probabilities"); legend(["Ground State", "1st Excited State", "2nd Excited State"]); grid()
```
![3-level NOT gate](img/3levelNOT.svg)

Instead of working perfectly, our system leaks into the 2nd energy level (we'll quantify this [later](#fidelity)). This becomes problematic when we're trying perform useful computations. DRAG is one possible remedy where introducing drive detuning and a second quadrature control, that is $$ℇ(t)=ℇ<sup>x</sup>(t)cos(ω<sub>d</sub>t)+ℇ<sup>y</sup>(t)sin(ω<sub>d</sub>t)$$, helps eliminate some of the leakage. Let's see how much of an improvement we get. We'll need to define a few more functions for the new controls.

```jldoctest example1
function DRAGx(t::Real,p::Vector)
    # Ɛˣ for fifth order DRAG
    tg  = p[1] # gate time
    σ   = p[2] # standard deviation
    A   = p[3] # amplitude
    Δ   = p[4] # anharmonicity
    λ   = p[5] # relative strength of transitions
    Ɛˣ = rotgaussianpulse(t,[tg,σ,A]) + (λ^2-4)*rotgaussianpulse(t,[tg,σ,A])^3/(8*Δ^2) - (13λ^4-76λ^2+112)*rotgaussianpulse(t,[tg,σ,A])^5/(128Δ^4)
    return Ɛˣ
end

function DRAGy(t::Real,p::Vector)
    # Ɛʸ for fifth order DRAG
    tg  = p[1] # gate time
    σ   = p[2] # standard deviation
    A   = p[3] # amplitude
    Δ   = p[4] # anharmonicity
    λ   = p[5] # relative strength of transitions
    B   = inv(√(2π)*σ*erf(tg/(√(8)*σ))-tg*gaussian(0.5tg,σ))
    Ɛˣ′ = A*B*Schrodinger.gaussianprime(t,σ)
    Ɛʸ = -Ɛˣ′/Δ + 33*(λ^2-2)*rotgaussianpulse(t,[tg,σ,A])^2*Ɛˣ′/(24*Δ^3)
    return Ɛʸ
end

function detuning(t::Real,p::Vector)
    # dynamical detuning for fifth order DRAG
    tg  = p[1] # gate time
    σ   = p[2] # standard deviation
    A   = p[3] # amplitude
    Δ   = p[4] # anharmonicity
    λ   = p[5] # relative strength of transitions
    δ₁ = (λ^2-4)*rotgaussianpulse(t,[tg,σ,A])^2/(4*Δ) - (λ^4-7λ^2+12)*rotgaussianpulse(t,[tg,σ,A])^4/(16Δ^3)
    return δ₁
end

Π₁ = basis(3,1) * basis(3,1)' # projector for 1st level
Hdet = Π₁ # affected by dynamical detuning
Hdx = create(3)/2 + destroy(3)/2 # affected by Ɛˣ
Hdy = im*create(3)/2 - im*destroy(3)/2 # affected by Ɛʸ

res3 = sesolve([Hc,(Hdet,detuning,[tg,σ,π,Δ,λ]),(Hdx,DRAGx,[tg,σ,π,Δ,λ]),(Hdy,DRAGy,[tg,σ,π,Δ,λ])],g,tspan)
figure()
plot(res3.times*1e9,levelprobs(res3.states)); xlabel("Time (ns)"); ylabel("Level Probabilities"); legend(["Ground State", "1st Excited State", "2nd Excited State"]); grid()
```
![3-level DRAG](img/3levelDRAG.svg)

There is a noticeable improvement over the simple Gaussian pulse, but how much better is the new gate?

## Fidelity

We can calculate the fidelity of our gates by comparing their output to the ideal case. Our gates behave ideally when $$λ=0$$ and there is no leakage into the 2nd excited state. Using the same measure of error as in \[[1]], we can take the overall gate fidelity to be the average of gate fidelities when using the 6 axial states on the Bloch sphere as inputs.

```jldoctest example1
tgs = (3:9)*1e-9 # gate times
Fg_res = Matrix{Float64}(length(tgs),2) # initialize matrix for number of solves
axialStates  = [normalize!(Ket([1,1,0])),   # +X
                normalize!(Ket([1,-1,0])),  # -X
                normalize!(Ket([1,im,0])),  # +Y
                normalize!(Ket([1,-im,0])), # -Y
                normalize!(Ket([1,0,0])),   # +Z
                normalize!(Ket([0,1,0]))]   # -Z
axialOperators = [] # need density operators too
for i = 1:6
    push!(axialOperators,axialStates[i]*axialStates[i]')
end
Uideal = complex(qzero(3)); Uideal[1:2,1:2] = data(σx)
for (i,tg) in enumerate(tgs)
    sum1 = 0 # sum of Gaussian gate fidelities
    sum2 = 0 # sum of DRAG gate fidelities
    for j = 1:6
        # Gaussian
        res4_1 = sesolve([Hc,(Hd,rotgaussianpulse,[tg,σ,π])],axialStates[j],(-tg/2,tg/2))
        sum1 += trace(Uideal*axialOperators[j]*Uideal'*(res4_1.states[end]*res4_1.states[end]'))
        # DRAG
        res4_2 = sesolve([Hc,(Hdet,detuning,[tg,σ,π,Δ,λ]),(Hdx,DRAGx,[tg,σ,π,Δ,λ]),(Hdy,DRAGy,[tg,σ,π,Δ,λ])],axialStates[j],(-tg/2,tg/2))
        sum2 += trace(Uideal*axialOperators[j]*Uideal'*(res4_2.states[end]*res4_2.states[end]'))
    end
    Fg_res[i,:] = [sum1/6 sum2/6] # take average
end
figure()
plot(tgs*1e9,1.-Fg_res); ylim([10e-8,1]); title("Average gate fidelity averaging over all input states"); yscale("log"); xlabel("Gate Time (ns)"); ylabel("Gate Error 1-Fg"); legend(["Gaussian","DRAG 5th Order"]); grid()
```
![NOT-gate fidelities](img/fidelities.svg)

For a gate time of 6ns, taking advantage of DRAG results in a gate error that is 2 orders of magnitude less than when using Gaussian pulses.

## References

\[[1]] F. Motzoi, J.M. Gambetta, P. Rebentrost, and F. K. Wilhelm, "Simple pulses for elimination of leakage in weakly nonlinear qubits," Phys. Rev. Lett. 103, 110501 (2009).

[1]: http://dx.doi.org/10.1103/PhysRevLett.103.110501
