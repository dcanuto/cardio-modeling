function modelhemo!(system::CVSystem,n::Int64,ID::Int64)
    system.branches.A[ID][n+2,system.solverparams.JL] = (system.hemo.Ph/
        system.branches.beta[ID][end]+sqrt(system.branches.A0[ID][end]))^2;
    system.branches.W2[ID] = system.branches.W1[ID]-
        8*((system.branches.A[ID][n+2,system.solverparams.JL]*(system.branches.beta[ID][end]/
        (2*system.solverparams.rho))^2)^0.25-system.branches.c0[ID][end]);
    system.branches.Q[ID][n+2,system.solverparams.JL] = 0.5*
        system.branches.A[ID][n+2,system.solverparams.JL]*(system.branches.W1[ID]+
        system.branches.W2[ID])

    system.hemo.Vloss += system.branches.Q[ID][n+2,system.solverparams.JL]*
        system.solverparams.h;

    for i = 1:length(system.branches.children[ID])
        cID = system.branches.children[ID][i];
        system.branches.A[cID][n+2,1] = (system.hemo.Ph/system.branches.beta[cID][end]+
            sqrt(system.branches.A0[cID][end]))^2;
        system.branches.W1[cID] = system.branches.W2[ID]+8*((system.branches.A[cID][n+2,1]*
            (system.branches.beta[cID][end]/(2*system.solverparams.rho))^2)^0.25-
            system.branches.c0[cID][end]);
        system.branches.Q[cID][n+2,1] = 0.5*system.branches.A[cID][n+2,1]*
            (system.branches.W1[cID]+system.branches.W2[cID]);
        system.hemo.Vloss += system.branches.Q[cID][n+2,1]*system.solverparams.h;
    end

end
