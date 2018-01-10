function fdist(x::Vector{Float64},system::CVSystem,n::Int64,ID::Int64)
    # non-dimensionalizing parameters
    Vs = system.branches.term[ID].V[n+1,1];
    vs = system.branches.c0[ID][end];
    ts = system.branches.k[ID]/vs;

    f1 = x[2]-system.branches.term[ID].V0[1]/Vs-system.branches.beta[ID][end]*
        system.branches.term[ID].C[1]/Vs*(vs^2*(2*system.solverparams.rho/
        system.branches.beta[ID][end])*((system.branches.W1end[ID]/vs-x[1])/8+
        system.branches.c0[ID][end]/vs)^2-sqrt(system.branches.A0[ID][end]));
    f2 = (x[2]-system.branches.term[ID].V[n+1,1]/Vs)/system.solverparams.h*ts-
        ts/Vs*(vs^5*(2*system.solverparams.rho/system.branches.beta[ID][end])^2*
        ((system.branches.W1end[ID]/vs-x[1])/8+system.branches.c0[ID][end]/vs)^4*
        0.5*(system.branches.W1end[ID]/vs+x[1]))+ts/Vs*system.branches.term[ID].Q[n+1,1];

    f = [f1,f2];
end
