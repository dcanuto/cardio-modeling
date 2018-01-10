function newtondist!(system::CVSystem,n::Int64,ID::Int64)
    # non-dimensionalizing parameters
    Vs = system.branches.term[ID].V[n+1,1];
    vs = system.branches.c0[ID][end];
    ts = system.branches.k[ID]/vs;

    # initial guess
    x0 = zeros(2);
    x0[2] = system.branches.term[ID].V[n+1,1];
    Pest = (1/system.branches.term[ID].C[1]*
        (x0[2] - system.branches.term[ID].V0[1]));
    Aest = (Pest/system.branches.beta[ID][end]+
        sqrt(system.branches.A0[ID][end]))^2;
    x0[1] = (system.branches.W1end[ID] + 8*(
        system.branches.c0[ID][end]-
        sqrt(system.branches.beta[ID][end]/(2*
        system.solverparams.rho))*Aest^0.25));

    # non-dimensionalize
    x0[1] = x0[1]/vs;
    x0[2] = x0[2]/Vs;

    # assign function/Jacobian handles
    f = fdist;
    J = Jdist;

    # setup for iterations
    xx = x0;
    x = zeros(length(x0));
    N = 1;

    # max step size for line searches
    stpmax = 100*max(sqrt(norm(xx)),length(x0));
    # println(stpmax)

    # Newton iterations
    while N <= system.solverparams.maxiter
        # determine Jacobian, check invertibility
        JJ = J(xx,system,n,ID);
        D = diagm(maximum!(zeros(length(xx)),abs(JJ)).^-1);
        # println(JJ)
        # println(D)
        # println(D*JJ)
        # println(cond(D*JJ))
        # println(det(D*JJ))
        if abs(det(D*JJ)) < system.solverparams.epsJ
            println(JJ)
            println(D*JJ)
            print(det(D*JJ))
            error("Distal Newton Jacobian is singular.");
        end
        # compute gradient of line search objective function
        fvec = D*f(xx,system,n,ID);
        # println(f(xx,system,n,state))
        # println(fvec)
        fold = 0.5*dot(fvec,fvec);
        g = (D*JJ)'*fvec;
        # compute newton step
        # println(inv(D*JJ))
        s = -inv(D*JJ)*fvec;
        # line search to update state vector
        fn, xn, check = linedist(xx,fold,g,s,stpmax,system,n,ID);
        # println(xn[1]*vs)
        # check if sufficiently close to root
        JJ = J(xn,system,n,ID);
        D = diagm(maximum!(zeros(length(xx)),abs(JJ)).^-1);
        fvec = D*f(xn,system,n,ID);
        # println(fvec)
        if norm(fvec) <= system.solverparams.epsN
            x = xn;
            x[1] = x[1]*vs;
            x[2] = x[2]*Vs;
            # println(x)
            # println(system.branches.W2root)
            # println(inv(D*JJ))
            # println(fvec)
            system.solverparams.totaliter+=N;
            break
        end
        if check
            test = 0.
            den = max(fn,0.5*length(xn));
            for i = 1:length(xn)
                temp = abs(g[i])*max(abs(xn[i]),1.)/den;
                if temp > test
                    test = temp;
                end
            end
            if test < 1e-6
                check = true;
                println("Warning: gradient of objective function close to zero.")
            else
                check = false;
            end
        end
        # check for divergence
        if norm(fvec,Inf) >= system.solverparams.maxval
            system.solverparams.totaliter+=N;
            println(xn)
            println(fvec)
            error("Newton iteration diverged.");
        end
        N+=1;
        xx = xn;
        if N == system.solverparams.maxiter
            println(JJ)
            # println(D)
            # println(D*JJ)
            println(xn)
            println(fvec)
            println(norm(fvec))
            error("Newton iteration failed to converge.");
        end
    end

    # update based on converged solution
    system.branches.term[ID].V[n+2,1] = x[2];
    system.branches.term[ID].P[n+2,1] = 1/system.branches.term[ID].C[1]*(x[2]-
        system.branches.term[ID].V0[1]);
    system.branches.A[ID][n+2,system.solverparams.JL] = (2*system.solverparams.rho/
        system.branches.beta[ID][end])^2*((system.branches.W1end[ID]-x[1])/8+
        system.branches.c0[ID][end])^4;
    system.branches.Q[ID][n+2,system.solverparams.JL] =
        system.branches.A[ID][n+2,system.solverparams.JL]*0.5*
        (system.branches.W1end[ID]+x[1]);
end
