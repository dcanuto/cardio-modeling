function coupleproximal!(system::CVSystem,n::Int64)
    j = Int8(1);
    omega = Float64(0.2);
    Qest = Float64;
    Vest = Float64;
    Pest = Float64;
    Pao = Float64;
    Aest = Float64;
    Uest = Float64;
    Qnew = Float64;
    W1 = Float64;
    iterr = Float64;

    if (system.heart.lv.P[n+1]/mmHgToPa > system.branches.P[1][n+1,1] ||
        system.branches.Q[1][n+1,1] > 0)
        AV = "open"
    else
        AV = "closed"
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
        # coupling loop
        while j <= system.solverparams.maxiter
            if j == 1
                Qest = system.branches.Q[1][n+1,1];
            end
            Vest = system.heart.lv.V[n+1] - system.solverparams.h*Qest;
            Pest = system.heart.lv.E[n+2]*(Vest - system.heart.lv.V0);
            Pao = Pest - Qest*system.branches.Rao;
            Aest = (Pao/system.branches.beta[1][end]+
                sqrt(system.branches.A0[1][end]))^2;
            W1 = (system.branches.W2root + 8*(sqrt(system.branches.beta[1][end]/
                (2*system.solverparams.rho))*Aest^0.25-
                system.branches.c0[1][end]));
            Uest = 0.5*(W1 + system.branches.W2root);
            Qnew = (1-omega)*Qest + omega*Aest*Uest;
            iterr = abs(Qnew - Qest);
            if iterr <= system.solverparams.abserror
                system.branches.A[1][n+2,1] = Aest;
                system.branches.Q[1][n+2,1] = Qest;
                system.heart.lv.V[n+2] = Vest;
                break
            else
                Qest = Qnew;
                j+=1
            end
            if j == system.solverparams.maxiter
                error("Proximal coupling failed to converge.")
            end
        end
    end
end
