function Jdist(x::Vector{Float64},system::CVSystem,n::Int64,ID::Int64)
    # non-dimensionalizing parameters
    Vs = system.branches.term[ID].V[n+1,1];
    vs = system.branches.c0[ID][end];
    ts = system.branches.k[ID]/vs;

    J11 = system.branches.term[ID].C[1]*system.solverparams.rho*vs^2/(2*Vs)*
        ((system.branches.W1end[ID]/vs-x[1])/8+system.branches.c0[ID][end]/vs)^2;
    J12 = 1;
    J21 = 2*ts*vs^5*system.solverparams.rho^2/(Vs*system.branches.beta[ID][end]^2)*
        (((system.branches.W1end[ID]/vs-x[1])/8+system.branches.c0[ID][end]/vs)^3*0.5*
        (system.branches.W1end[ID]/vs+x[1])-((system.branches.W1end[ID]/vs-x[1])/8+
        system.branches.c0[ID][end]/vs)^4)
    J22 = ts/system.solverparams.h;
    J1 = [J11 J12];
    J2 = [J21 J22];

    J = [J1;J2];
end
