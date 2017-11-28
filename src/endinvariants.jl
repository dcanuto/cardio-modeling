function endinvariants!(system::CVSystem,n::Int64)
    lf = Float64;
    W1 = Vector{Float64}[];
    dW1 = Float64;
    dQ = Float64;
    dA = Float64;
    Qint = Float64;
    Aint = Float64;

    for i = 1:length(system.branches.ID)
        if isempty(system.branches.children[i])
            # forward wave speed
            lf = (system.branches.Q[i][n+1,system.solverparams.JL]/
                system.branches.A[i][n+1,system.solverparams.JL] +
                sqrt(0.5*system.branches.beta[i][end]/system.solverparams.rho)*
                system.branches.A[i][n+1,system.solverparams.JL]^0.25);
            # right-running invariants at current time
            W1 = ([system.branches.Q[i][n+1,system.solverparams.JL-1]/
                system.branches.A[i][n+1,system.solverparams.JL-1]+
                4*(sqrt(0.5*system.branches.beta[i][end]/
                system.solverparams.rho)*
                system.branches.A[i][n+1,system.solverparams.JL-1]^0.25 -
                system.branches.c0[i][end])]);
            push!(W1,system.branches.Q[i][n+1,system.solverparams.JL]/
                system.branches.A[i][n+1,system.solverparams.JL]+
                4*(sqrt(0.5*system.branches.beta[i][end]/
                system.solverparams.rho)*
                system.branches.A[i][n+1,system.solverparams.JL]^0.25 -
                system.branches.c0[i][end]));
            dW1 = (W1[2] - W1[1])/system.branches.k[i];
            # change in solution variables
            dQ = ((system.branches.Q[i][n+1,system.solverparams.JL-1] -
                system.branches.Q[i][n+1,system.solverparams.JL])/
                system.branches.k[i]);
            dA = ((system.branches.A[i][n+1,system.solverparams.JL-1] -
                system.branches.A[i][n+1,system.solverparams.JL])/
                system.branches.k[i]);
            # interpolated solution variables
            Qint = (system.branches.Q[i][n+1,system.solverparams.JL-1] -
                dQ*lf*system.solverparams.h);
            Aint = (system.branches.A[i][n+1,system.solverparams.JL-1] -
                dA*lf*system.solverparams.h);
            # update right-running invariant by linear interp.
            system.branches.W1end[i] = (W1[2] - dW1*lf*system.solverparams.h -
                system.solverparams.diffusioncoeff*system.solverparams.mu/
                system.solverparams.rho*system.solverparams.h*(Qint/Aint^2));
        end
    end
end
