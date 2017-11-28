function fdouble(x::Vector{Float64},system::CVSystem,parentID::Int64,children::Vector{Int64})
    # conservation of mass
    f1 = x[1] - x[3] - x[5];

    # total P
    P1 = system.branches.beta[parentID][end]*(sqrt(x[2]) - sqrt(
        system.branches.A0[parentID][end]));
    P2 = system.branches.beta[children[1]][end]*(sqrt(x[4]) - sqrt(
        system.branches.A0[children[1]][end]));
    P3 = system.branches.beta[children[2]][end]*(sqrt(x[6]) - sqrt(
        system.branches.A0[children[2]][end]));

    f2 = (P1 + 0.5*system.solverparams.rho*(x[1]/x[2])^2)-
        (P2 + 0.5*system.solverparams.rho*(x[3]/x[4])^2);
    f3 = (P1 + 0.5*system.solverparams.rho*(x[1]/x[2])^2)-
        (P3 + 0.5*system.solverparams.rho*(x[5]/x[6])^2);

    # Riemann invariants
    f4 = (x[1]/x[2]) + 4*(sqrt(0.5*system.branches.beta[parentID][end]/
        system.solverparams.rho)*x[2]^0.25 - system.branches.c0[parentID][end])-
        system.branches.W1[parentID];
    f5 = (x[3]/x[4]) - 4*(sqrt(0.5*system.branches.beta[children[1]][end]/
        system.solverparams.rho)*x[4]^0.25 -
        system.branches.c0[children[1]][end]) - system.branches.W2[children[1]];
    f6 = (x[5]/x[6]) - 4*(sqrt(0.5*system.branches.beta[children[2]][end]/
        system.solverparams.rho)*x[6]^0.25 -
        system.branches.c0[children[2]][end]) - system.branches.W2[children[2]];

    f = [f1,f2,f3,f4,f5,f6];
end
