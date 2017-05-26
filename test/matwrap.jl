using BenchmarkTools, Schrodinger
n = 10
A = operator(rand(n,n),(n,),true);
B = operator(rand(n,n),(n,),true);
A*B;
A.data*B.data;
@time for i=1:300000;A*B;end
@time for i=1:300000;A.data*B.data;end

function test_operator()
    for n in [10,16,30]
        A = operator(rand(n,n),(n,),true)
        B = operator(rand(n,n),(n,),true)
        @show @benchmark $A*$B
        @show @benchmark $A.data*$B.data
    end
end

test_operator()

@benchmark begin
N = 10
omega_a = 1.0
omega_c = 1.0
g = 0.05
a = tensor(sigma0, destroy(N))
sm = tensor(sigmaplus, qeye(N))
sz = tensor(sigmaz, qeye(N))
smdsm = (sm'*sm)
H = 0.5 * omega_a * sz + omega_c * a' * a + g * (a' * sm + a * sm')
psi0 = tensor(basis(2,1),basis(N,0))
end

#python -m timeit -s "import qutip as qt" "N = 10; omega_a = 1.0; omega_c = 1.25; g = 0.05; a = qt.tensor(qt.identity(2), qt.destroy(N)); sm = qt.tensor(qt.destroy(2), qt.identity(N)); sz = qt.tensor(qt.sigmaz(), qt.identity(N)); H = 0.5 * omega_a * sz + omega_c* a.dag() * a + g *(a.dag() * sm + a * sm.dag())"
