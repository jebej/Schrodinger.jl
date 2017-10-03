using Schrodinger, PyPlot

function opt_pihalfx(n)
    Hd = qzero(2)
    Hc = [π*σx, π*σy]
    t = 1
    #u_init = rand(n,1) .- 0.5
    u_init = zeros(n,2)
    Ut = expm(-1im*π*σx/4)
    # Create objective function type
    O = NormPSU(Ut,Hd,Hc,t,n)
    grape(O,u_init)
end

function opt_hadamard(n)
    Hd = σz
    Hc = [σx]
    t = 10.0
    #u_init = rand(n,1) .- 0.5
    u_init = zeros(n,1)
    Ut = Operator([1/√2 1/√2; 1/√2 -1/√2])
    # Create objective function type
    O = NormPSU(Ut,Hd,Hc,t,n)
    grape(O,u_init)
end

function opt_3lvlNOT(n)
    Δ = 2π*(-400) # anharmonicity
    Hd = Δ*Operator(basis(3,2)) # drift Hamiltonian
    Hc = [create(3)/2+destroy(3)/2, im*create(3)/2-im*destroy(3)/2]
    t = 3
    u_init = [-Δ*Schrodinger.gaussianpulse.(linspace(-t/2,t/2,n),[[t/2,t,0,0,pi]]) linspace(-Δ/2,Δ/2,n)]
    Ut = Operator([0 1 0; 1 0 0; 0 0 1]) # 3lvl NOT gate
    # Create objective function type
    O = CoherentSubspaces(Ut,1:2,Hd,Hc,t,n) # care only about computational subspace
    grape(O,u_init)
end

function opt_2Q_QFT(n)
    Hd = 0.5 * (σx⊗σx + σy⊗σy + σz⊗σz)
    Hc = [0.5*σx⊗σ0, 0.5*σy⊗σ0, 0.5*σ0⊗σx, 0.5*σ0⊗σy]
    t = 6.0
    u_init = zeros(n,4)
    #u_init = rand(n,4) .- 0.5
    Ut = Operator([0.5  0.5    0.5  0.5
                   0.5  0.5im -0.5 -0.5im
                   0.5 -0.5    0.5 -0.5
                   0.5 -0.5im -0.5  0.5im],
                  (2,2))
    # Create objective function type
    O = NormPSU(Ut,Hd,Hc,t,n)
    grape(O,u_init)
end

function plotgrape(tf,uf)
    figure()
    step(linspace(0,tf,size(uf,1)+1),[uf[1:1,:]; uf])
    legend(["Control $i" for i in 1:size(uf,2)])
    tight_layout()
    grid()
end

benchgrape(O,u,n=100) = for i=1:n; grape(O,u); end
