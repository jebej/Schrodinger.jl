```@meta
DocTestSetup  = quote
    using Schrodinger, PyPlot
end
```

```@setup driven_TLS
!isdir("img") && mkdir("img")
```

# The Driven Two Level System

The driven two level system (TLS) is one of the very few non-trivial time dependent problem whose dynamics can be solved for analytically. As such it provides a good benchmark for numerical solvers, since we can easily quantify their performance.

## The Hamiltonian

The Hamiltonian of a two level system is usually expressed with a $$\hat{σ}_z$$ Pauli operator, but for simplicity, we use a normal number operator, which we multiply by the frequency of the transition between the ground state and the excited state: $$\hat{H}_0 = \hbar ω \hat{n}$$.

The drive consists of an oscillating field in a direction transverse to the TLS. We use a circularly polarized field, meaning that the component of the drive for the $$0 \rightarrow 1$$ transition is complex and is the complex conjugate of the component of the drive for the $$1 \leftarrow 0$$ transition. Explicitly, if the drive has strength $$\gamma$$ (in angular frequency), the drive matrix is:

\begin{align\*}
\hat{H}_d = \hbar \gamma \begin{pmatrix}
0 & e^{i \omega_d t} \\\\
e^{-i \omega_d t} & 0
\end{pmatrix},
\end{align\*}
where $$\omega_d$$ is the drive frequency.

Note that if the driving field is not circularly polarized, for example if the drive is $$\hbar \gamma \hat{\sigma}_x \cos(\omega_d t)$$ (corresponding to linear polarization), the dynamics of the system are no longer solvable analytically.

## Numerics

Before starting the mathematical derivation, let's let the computer solve the problem for us:

```@example driven_TLS
using Schrodinger, PyPlot

# frequencies, starting state and number operator
ω = 1*2π
ωd = 0.9*2π
γ = π/2
g = basis(2,0)
n = numberop(2)

# drive functions, two for circular polarization (cis(x) == exp(1im*x))
# p[1] is the frequency, and p[2] the drive strength in angular frequency
f1m = (t,p) -> p[2]*cis(p[1]*t)
f1p = (t,p) -> p[2]*cis(-p[1]*t)

# solve! we use a drive strength of γ
r1 = sesolve((ω*n, (σ₋,f1m,[ωd,γ]), (σ₊,f1p,[ωd,γ])), g, (0.0,3.0), saveat=0.01)

plot(r1.times,levelprobs(r1.states)); grid();
xlabel("Time"); ylabel("Level Probabilities");
legend(["|g⟩", "|e⟩"]);
tight_layout(true); savefig(joinpath("img","driven_TLS_1.svg")); # hide
```
![driven_TLS_1](img/driven_TLS_1.svg)

Notice that we chose a drive frequency below that of the TLS frequency. This results in a TLS oscillation frequency *higher* than π (the drive strength), however, the e never gets fully occupied.

## Derivation

We now derive the dynamics of the system analytically, so that we may verify our numerics.

First, we move to the interaction picture in order to ignore the precession about the z-axis. We are thus left only with $$\hat{H}_I = e^{i \hat{H}_0 t / \hbar} \hat{H}_D e^{-i \hat{H}_0 t / \hbar}$$. For a state with time dependent coefficients $$c_0$$ and $$c_1$, the Schrödinger equation thus reads:

\begin{align\*}
i \hbar \dot{|\psi\rangle}_I &= \hat{H}_I |\psi\rangle_I, \\\\
i \hbar \begin{pmatrix} \dot{c}_0 \\\\ \dot{c}_1 \end{pmatrix} &= \hbar \gamma \begin{pmatrix} 0 & e^{i \Delta\omega t} \\\\ e^{-i \Delta\omega t} & 0 \end{pmatrix} \begin{pmatrix} c_0 \\\\ c_1 \end{pmatrix}, \\\\
 &= \hbar \gamma \begin{pmatrix} e^{i \Delta\omega t} c_1 \\\\ e^{-i \Delta\omega t} c_0 \end{pmatrix},
\end{align\*}
where $$\Delta\omega = \omega_d - \omega$$.

We thus have a system of coupled differential equations to solve, with initial condition $$c_0(t=0) = 1$$ and normalization $$|c_0|^2 + |c_1|^2 = 1$$. We first decouple the system by isolating $$c_0(t)$$ in the second equation, differentiating it and inserting the result in the second equation:

\begin{align\*}
c_0(t) &= i/\gamma  e^{i \Delta\omega t} \dot{c}_1, \\\\
\dot{c}_0(t) &= i/\gamma  \left( e^{i \Delta\omega t} \ddot{c}_1 + i \Delta\omega e^{i \Delta\omega t} \dot{c}_1 \right),
\end{align\*}

\begin{align\*}
-1/\gamma  \left( e^{i \Delta\omega t} \ddot{c}_1 + i \Delta\omega e^{i \Delta\omega t} \dot{c}_1 \right) &= \gamma  e^{i \Delta\omega t} c_1, \\\\
e^{i \Delta\omega t} \ddot{c}_1 + i \Delta\omega e^{i \Delta\omega t} \dot{c}_1 &= -\gamma^2  e^{i \Delta\omega t} c_1.
\end{align\*}

Finally, since the equation must be valid for all $$t$$, we drop the exponential to obtain a homogeneous second order differential equation:

\begin{align\*}
\ddot{c}_1 + i \Delta\omega \dot{c}_1 + \gamma^2 c_1 = 0.
\end{align\*}

At this point we guess a solution $$ c_1(t) \propto e^{i \alpha t}$$ which, when inserted in the differential equation, gives a 2nd order polynomial in $$\alpha$$:

\begin{align\*}
-\alpha^2 c_1 + i \Delta\omega i \alpha c_1 + \gamma^2 c_1 &= 0, \\\\
\alpha^2 + \Delta\omega \alpha - \gamma^2 &= 0,
\end{align\*}

with solutions

\begin{align\*}
\alpha_\pm &= - \Delta\omega/2 \pm \sqrt{ ( \Delta\omega/2 )^2 + \gamma^2 }, \\\\
           &= - \Delta\omega/2 \pm \Omega,
\end{align\*}
where we have introduced the well know *Rabi frequency* $$\Omega = \sqrt{ ( \Delta\omega/2 )^2 + \gamma^2 }$$.


As expected, a second order differential equation gives us two solutions, with which we form the general solution:

\begin{align\*}
c_1(t) = a e^{-i \Delta\omega t /2} e^{i \Omega t} + b e^{-i \Delta\omega t /2} e^{-i \Omega t}.
\end{align\*}

We can eliminate $$b$$ by using the initial condition $$c_1(t=0) = 0$$, giving us $$b = -a$$. This allows us to simplify the general solution:

\begin{align\*}
c_1(t) = a e^{-i \Delta\omega t /2} 2i \sin(\Omega t) = c e^{-i \Delta\omega t /2} \sin(\Omega t).
\end{align\*}

Finally, to determine the value of $$a$$, we insert the general solution back in the original differential equation for $$c_1(t)$$ and make use of the normalization property:

\begin{align\*}
\dot{c_1} &= c ( -i \Delta\omega/2 e^{-i \Delta\omega t /2} \sin(\Omega t) + \Omega e^{-i \Delta\omega t /2} \cos(\Omega t) ), \\\\
          &= \gamma e^{-i \Delta\omega t} c_0
\end{align\*}

equating the equations and taking the absolute square:

\begin{align\*}
e^{-i \Delta\omega t /2} c ( -i \Delta\omega/2 \sin(\Omega t) + \Omega \cos(\Omega t) ) &= \gamma e^{-i \Delta\omega t} c_0, \\\\
|c|^2 ( (\Delta\omega/2)^2 \sin^2(\Omega t) + \Omega^2 \cos^2(\Omega t) ) &= \gamma^2 |c_0|^2 = \gamma^2 ( 1 - |c_1|^2 ), \\\\
|c|^2 ( (\Delta\omega/2)^2 \sin^2(\Omega t) + \Omega^2 \cos^2(\Omega t) + \gamma^2 \sin^2(\Omega t) )&= \gamma^2, \\\\
|c|^2 ( \Omega^2 \sin^2(\Omega t) + \Omega^2 \cos^2(\Omega t) )&= \gamma^2, \\\\
|c|^2 \Omega^2 &= \gamma^2, \\\\
|c|^2 &= \gamma^2 / \Omega^2.
\end{align\*}

The final solution for $$c_1(t)$$ then, up to a constant global phase is:

\begin{align\*}
c_1(t) = \gamma / \Omega e^{-i \Delta\omega t /2} \sin(\Omega t).
\end{align\*}

This equation is quite interesting, and it perfectly explains the non-trivial behavior of the driven TLS. It shows the switching frequency, which is the much celebrated Rabi frequency, and also shows that the maximum transfer probability from the ground to excited state $$|c|^2$$ is a Lorentzian with respect to $$\Delta\omega$$, the detuning between the drive and the TLS frequency.

## Comparison

In order to compare to the simulation done at the beginning, we need to go back in the Schrödinger picture. This is done by multiplying the solution above by $$e^{-i ω t}$$. Alternatively, we could perform the simulation in the interaction picture. Now, let's plot the time evolution!

```@example driven_TLS
c₁(t,γ,Δω,ω) = γ/√(γ^2+(Δω/2)^2) * cis(-t*(Δω/2+ω)) * sin(√(γ^2+(Δω/2)^2)*t)

Δω = ωd-ω
c1 = c₁.(r1.times,γ,Δω,ω) # compute data for c₁
analytic = [1-abs2.(c1) abs2.(c1)] # 1-|c₁|² and |c₁|²

plot(r1.times,analytic)
legend(["|g⟩", "|e⟩", "|g⟩ (analytic)", "|e⟩ (analytic)"]);
savefig(joinpath("img","driven_TLS_2.svg")); # hide
```
![driven_TLS_2](img/driven_TLS_2.svg)

It worked! It is hard to tell that the two lines are on top of each other, but we can numerically compare the data:

```@example driven_TLS
analytic ≈ levelprobs(r1.states)
```

Note that computationally, it is much more advantageous to perform the simulation in the interaction picture since we avoid simulating the constant precession of the TLS around the z-axis.

```@example driven_TLS
# in the interaction picture
r2 = sesolve((qzero(2), (σ₋,f1m,[Δω,γ]), (σ₊,f1p,[Δω,γ])), g, (0.0,3.0), saveat=0.01)

analytic ≈ levelprobs(r2.states)
```

## Chevron Pattern

A cool experiment that can be done with qubits is to verify that the TLS dynamics hold over a range of driving frequencies. We drive a qubit at various frequencies and measure the excited state probability over time. This leads to a very nice "chevron pattern" plot:

```@example driven_TLS
t = linspace(0,12,101)
ωd_ratio = linspace(0.6,1.4,101).'
chev = abs2.(c₁.(t,γ,ω.*ωd_ratio-ω,ω))
figure()
imshow(chev,origin="lower",extent=[0.6,1.4,0,12],aspect=1/15)
xlabel("\$ω_d/ω\$"); ylabel("Time");
tight_layout(true); savefig(joinpath("img","driven_TLS_3.svg")); # hide
```
![driven_TLS_2](img/driven_TLS_3.svg)
