function rootinvariant!(system::CVSystem,n::Int64)
    lb = Float64;
    W2 = Vector{Float64}[];
    dW2 = Float64;
    dQ = Float64;
    dA = Float64;
    Qint = Float64;
    Aint = Float64;

    # backward wave speed at proximal end
    lb = (0.5*((system.branches.Q[1][n+1,2]/system.branches.A[1][n+1,2]-
        sqrt(0.5*system.branches.beta[1][end]/system.solverparams.rho)*
        system.branches.A[1][n+1,2]^0.25)+
        (system.branches.Q[1][n+1,1]/system.branches.A[1][n+1,1]-
        sqrt(0.5*system.branches.beta[1][end]/system.solverparams.rho)*
        system.branches.A[1][n+1,1]^0.25)));
    # left-running invariant at current time
    W2 = ([system.branches.Q[1][n+1,1]/system.branches.A[1][n+1,1]-
        4*(sqrt(0.5*system.branches.beta[1][end]/system.solverparams.rho)*
        system.branches.A[1][n+1,1]^0.25 - system.branches.c0[1][end])]);
    push!(W2,system.branches.Q[1][n+1,2]/system.branches.A[1][n+1,2]-
        4*(sqrt(0.5*system.branches.beta[1][end]/system.solverparams.rho)*
        system.branches.A[1][n+1,2]^0.25 - system.branches.c0[1][end]))
    dW2 = (W2[2] - W2[1])/system.branches.k[1];
    # change in solution variables
    dQ = ((system.branches.Q[1][n+1,2] - system.branches.Q[1][n+1,1])/
        system.branches.k[1]);
    dA = ((system.branches.A[1][n+1,2] - system.branches.A[1][n+1,1])/
        system.branches.k[1]);
    # interpolated state variables
    Qint = system.branches.Q[1][n+1,1] - dQ*lb*system.solverparams.h;
    Aint = system.branches.A[1][n+1,1] - dA*lb*system.solverparams.h;
    # advance proximal boundary invariant
    system.branches.W2root = (W2[1] - dW2*lb*system.solverparams.h-
        system.solverparams.diffusioncoeff*system.solverparams.mu/
        system.solverparams.rho*system.solverparams.h*(Qint/Aint^2));
end
