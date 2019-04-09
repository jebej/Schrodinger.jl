using Schrodinger, PyPlot
using Schrodinger: gaussianpulse

function opt_pihalfx(n)
    Hd = qzero(2)
    Hc = [π*σx, π*σy]
    t = 1
    #ui = rand(n,1) .- 0.5
    ui = zeros(n,2)
    Ut = exp(-1im*π*σx/4)
    # Create objective function type
    O = NormPSU(Ut,Hd,Hc,t,n)
    grape(O,ui)
end

function opt_hadamard(n)
    Hd = σz
    Hc = [σx]
    t = 10.0
    #ui = rand(n,1) .- 0.5
    ui = zeros(n,1)
    Ut = Operator([1/√2 1/√2; 1/√2 -1/√2])
    # Create objective function type
    O = NormPSU(Ut,Hd,Hc,t,n)
    grape(O,ui)
end

function opt_RF_3lvlNOT(n)
    Δ = 2π*(-0.4) # anharmonicity (GHz)
    Hd = Δ*Operator(basis(3,2)) # drift Hamiltonian
    Hc = [create(3)/2+destroy(3)/2, im*create(3)/2-im*destroy(3)/2]
    t = 3 # (ns)
    tvec = LinRange(-t/2,t/2,n)
    ui = [-Δ*gaussianpulse.(tvec,[[t/2,t,0,0,pi]]) LinRange(-Δ/2,Δ/2,n)]
    Ut = Operator([0 1 0; 1 0 0; 0 0 1]) # 3lvl NOT gate
    # Create objective function type
    O = CoherentSubspaces(Ut,1:2,Hd,Hc,t,n) # care only about computational subspace
    grape(O,ui)
end

function opt_3lvlNOT(n)
    ω = 2π*7.0 # 0-1 transition freq (GHz)
    Δ = 2π*(-0.4) # anharmonicity
    Hd = ω*Operator(basis(3,1))+(2ω+Δ)*Operator(basis(3,2)) # drift Hamiltonian
    Hc = [create(3)+destroy(3)]
    t = 3 # (ns)
    tvec = LinRange(-t/2,t/2,n)
    ui = gaussianpulse.(tvec,[[t/2,t,ω,0,pi]])
    #ui = 2 .* cos.(ω.*tvec)
    #Ut = Operator([0 1 0; 1 0 0; 0 0 1]) # 3lvl NOT gate
    Ut = (Operator([0 1 0; 1 0 0; 0 0 0]),Operator([0 0 0; 0 0 0; 0 0 -1]))
    # Create objective function type
    O = CoherentSubspaces(Ut,[1],Hd,Hc,t,n) # care only about computational subspace
    grape(O,ui)
end

function opt_2Q_QFT(n)
    Hd = 0.5 * (σx⊗σx + σy⊗σy + σz⊗σz)
    Hc = [0.5*σx⊗σ0, 0.5*σy⊗σ0, 0.5*σ0⊗σx, 0.5*σ0⊗σy]
    t = 6.0
    #ui = zeros(n,4)
    #ui = rand(n,4) .- 0.5
    ui = Schrodinger.randomsmooth(n,4)
    Ut = Operator([0.5  0.5    0.5  0.5
                   0.5  0.5im -0.5 -0.5im
                   0.5 -0.5    0.5 -0.5
                   0.5 -0.5im -0.5  0.5im],
                  (2,2))
    # Create objective function type
    O = NormPSU(Ut,Hd,Hc,t,n)
    grape(O,ui)
end

benchgrape(O,u,n=100) = for i=1:n; grape(O,u); end

function run_examples()
    r1 = opt_pihalfx(20)
    r2 = opt_hadamard(200)
    r3 = opt_RF_3lvlNOT(300)
    r4 = opt_3lvlNOT(1000)
    r5 = opt_2Q_QFT(500)
    rr = [r1, r2, r3, r4, r5]
    [plotgrape(r) for r in rr]
    return rr
end

function plotgrape(res::Schrodinger.GrapeResult)
    ui = res.ui
    uf = res.uf
    t = 0 : res.t/size(uf,1) : res.t
    figure()
    step(t,[uf[1:1,:]; uf])
    gca()[:set_prop_cycle](nothing)
    step(t,[ui[1:1,:]; ui],linestyle="dashed")
    legend(["Control $i" for i in 1:size(uf,2)])
    tight_layout(true)
    grid(true)
end
