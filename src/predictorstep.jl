function predictorstep!(system::CVSystem,n::Int64)
    for i = 1:length(system.branches.ID)
        system.branches.Aforward[i][:] = (0.5*
            (system.branches.A[i][n+1,system.solverparams.colsint+1] +
            system.branches.A[i][n+1,system.solverparams.colsint]) -
            system.solverparams.h/(2*system.branches.k[i])*
            (system.branches.Fp[i][system.solverparams.acolspre+1] -
            system.branches.Fp[i][system.solverparams.acolspre]));
        system.branches.Qforward[i][:] = (0.5*
            (system.branches.Q[i][n+1,system.solverparams.colsint+1] +
            system.branches.Q[i][n+1,system.solverparams.colsint]) -
            system.solverparams.h/(2*system.branches.k[i])*
            (system.branches.Fp[i][system.solverparams.qcolspre+1] -
            system.branches.Fp[i][system.solverparams.qcolspre]) -
            0.5*system.solverparams.h*system.solverparams.diffusioncoeff*π*
            system.solverparams.mu/system.solverparams.rho*
            (0.5*(system.branches.Q[i][n+1,system.solverparams.colsint+1]./
            system.branches.A[i][n+1,system.solverparams.colsint+1] +
            system.branches.Q[i][n+1,system.solverparams.colsint]./
            system.branches.A[i][n+1,system.solverparams.colsint])));
        system.branches.Abackward[i][:] = (0.5*
            (system.branches.A[i][n+1,system.solverparams.colsint] +
            system.branches.A[i][n+1,system.solverparams.colsint-1]) -
            system.solverparams.h/(2*system.branches.k[i])*
            (system.branches.Fp[i][system.solverparams.acolspre] -
            system.branches.Fp[i][system.solverparams.acolspre-1]));
        system.branches.Qbackward[i][:] = (0.5*
            (system.branches.Q[i][n+1,system.solverparams.colsint] +
            system.branches.Q[i][n+1,system.solverparams.colsint-1]) -
            system.solverparams.h/(2*system.branches.k[i])*
            (system.branches.Fp[i][system.solverparams.qcolspre] -
            system.branches.Fp[i][system.solverparams.qcolspre-1]) -
            0.5*system.solverparams.h*system.solverparams.diffusioncoeff*π*
            system.solverparams.mu/system.solverparams.rho*
            (0.5*(system.branches.Q[i][n+1,system.solverparams.colsint]./
            system.branches.A[i][n+1,system.solverparams.colsint] +
            system.branches.Q[i][n+1,system.solverparams.colsint-1]./
            system.branches.A[i][n+1,system.solverparams.colsint-1])));
    end
end
