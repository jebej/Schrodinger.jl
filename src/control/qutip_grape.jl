using PyCall
@pyimport qutip as qt
@pyimport qutip.qip.algorithms.qft as qft
@pyimport qutip.control.pulseoptim as cpo


function grape_qutip(Ut,Hd,Hc,N,t,n)
    # start point for the gate evolution
    U0 = qt.identity(N)
    res = cpo.optimize_pulse_unitary(Hd,Hc,U0,Ut,n,t,init_pulse_type="ZERO")
    uf = res[:final_amps]
    Uf = res[:evo_full_final]
    return uf, Uf, res
end

function opt_hadamard_qutip(n)
    Hd = qt.sigmaz()
    Hc = [qt.sigmax()]
    t = 10.0
    Ut = qt.hadamard_transform(1)
    grape_qutip(Ut,Hd,Hc,2,t,n)
end

function opt_2Q_QFT_qutip(n)
    Sx = qt.sigmax()
    Sy = qt.sigmay()
    Sz = qt.sigmaz()
    Si = qt.identity(2)
    # Drift Hamiltonian and the (four) control Hamiltonians
    H_d = (qt.tensor(Sx, Sx) + qt.tensor(Sy, Sy) + qt.tensor(Sz, Sz))
    H_c = [qt.tensor(Sx, Si), qt.tensor(Sy, Si), qt.tensor(Si, Sx), qt.tensor(Si, Sy)]
    # Target for the gate evolution - Quantum Fourier Transform gate
    U_targ = qft.qft(2)
    # Evolution time
    t = 6.0
    grape_qutip(U_targ,H_d,H_c,4,t,n)
end

#"import qutip as qt; import qutip.qip.algorithms.qft as qft; import qutip.control.pulseoptim as cpo; Sx = qt.sigmax();Sy = qt.sigmay();Sz = qt.sigmaz();Si = qt.identity(2);H_d = (qt.tensor(Sx, Sx) + qt.tensor(Sy, Sy) + qt.tensor(Sz, Sz));H_c = [qt.tensor(Sx, Si), qt.tensor(Sy, Si), qt.tensor(Si, Sx), qt.tensor(Si, Sy)];U_0 = qt.identity(4);U_targ = qft.qft(2);" "cpo.optimize_pulse_unitary(H_d, H_c, U_0, U_targ, 100, 6, init_pulse_type='ZERO')"
