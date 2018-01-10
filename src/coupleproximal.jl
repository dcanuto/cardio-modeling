function coupleproximal!(system::CVSystem,n::Int64)
    if (system.heart.lv.P[n+1]/mmHgToPa < system.branches.P[1][n+1,1] &&
        system.heart.av.zeta[n+1] == 0)
        AV = "closed";
    elseif (system.heart.av.zeta[n+1] >= 0 && system.heart.lv.P[n+1]/mmHgToPa >
        system.branches.P[1][n+1,1])
        AV = "opening";
    elseif system.heart.lv.P[n+1]/mmHgToPa < system.branches.P[1][n+1,1]
        AV = "closing";
    end

    if AV == "closed"
        # update left-running invariant by extrapolation
        rootinvariant!(system,n);
        # update right-running invariant
        system.branches.W1root = -system.branches.W2root;
        # update proximal A, Q w/ invariants
        system.branches.A[1][n+2,1] = ((2*system.solverparams.rho/
            system.branches.beta[1][end])^2*(0.125*(system.branches.W1root-
            system.branches.W2root) + system.branches.c0[1][end])^4);
        system.branches.Q[1][n+2,1] = (system.branches.A[1][n+2,1]*0.5*(
            system.branches.W1root + system.branches.W2root));
        system.heart.lv.V[n+2] = (system.heart.lv.V[n+1]+
            system.solverparams.h*(system.heart.la.Q[n+1]-
            system.branches.Q[1][n+1,1]));
    else
        # ventricular elastance at next time step
        elastancefn!(system,n+1);
        # update W2 by extrapolation
        rootinvariant!(system,n);
        # Newton iterations
        newtonav!(system,n,AV);
    end
end
