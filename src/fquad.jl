function fquad(x::Vector{Float64},system::CVSystem,parentID::Int64,children::Vector{Int64})
    # conservation of mass
    f1 = x[1] - x[3] - x[5] - x[7] - x[9];

    # total P
    P1 = system.branches.beta[parentID][end]*(sqrt(x[2]) - sqrt(
        system.branches.A0[parentID][end]));
    P2 = system.branches.beta[children[1]][end]*(sqrt(x[4]) - sqrt(
        system.branches.A0[children[1]][end]));
    P3 = system.branches.beta[children[2]][end]*(sqrt(x[6]) - sqrt(
        system.branches.A0[children[2]][end]));
    P4 = system.branches.beta[children[3]][end]*(sqrt(x[8]) - sqrt(
        system.branches.A0[children[3]][end]));
    P5 = system.branches.beta[children[4]][end]*(sqrt(x[10]) - sqrt(
        system.branches.A0[children[4]][end]));

    f2 = (P1 + 0.5*system.solverparams.rho*(x[1]/x[2])^2)-
        (P2 + 0.5*system.solverparams.rho*(x[3]/x[4])^2);
    f3 = (P1 + 0.5*system.solverparams.rho*(x[1]/x[2])^2)-
        (P3 + 0.5*system.solverparams.rho*(x[5]/x[6])^2);
    f4 = (P1 + 0.5*system.solverparams.rho*(x[1]/x[2])^2)-
        (P4 + 0.5*system.solverparams.rho*(x[7]/x[8])^2);
    f5 = (P1 + 0.5*system.solverparams.rho*(x[1]/x[2])^2)-
        (P5 + 0.5*system.solverparams.rho*(x[9]/x[10])^2);

    # Riemann invariants
    f6 = (x[1]/x[2]) + 4*(sqrt(0.5*system.branches.beta[parentID][end]/
        system.solverparams.rho)*x[2]^0.25 - system.branches.c0[parentID][end])-
        system.branches.W1[parentID];
    f7 = (x[3]/x[4]) - 4*(sqrt(0.5*system.branches.beta[children[1]][end]/
        system.solverparams.rho)*x[4]^0.25 -
        system.branches.c0[children[1]][end]) - system.branches.W2[children[1]];
    f8 = (x[5]/x[6]) - 4*(sqrt(0.5*system.branches.beta[children[2]][end]/
        system.solverparams.rho)*x[6]^0.25 -
        system.branches.c0[children[2]][end]) - system.branches.W2[children[2]];
    f9 = (x[7]/x[8]) - 4*(sqrt(0.5*system.branches.beta[children[3]][end]/
        system.solverparams.rho)*x[8]^0.25 -
        system.branches.c0[children[3]][end]) - system.branches.W2[children[3]];
    f10 = (x[9]/x[10]) - 4*(sqrt(0.5*system.branches.beta[children[4]][end]/
        system.solverparams.rho)*x[10]^0.25 -
        system.branches.c0[children[4]][end]) - system.branches.W2[children[4]];


    f = [f1,f2,f3,f4,f5,f6,f7,f8,f9,f10];
end
