function newtonav!(system::CVSystem,n::Int64,state::String)
    # non-dimensionalizing parameters
    Vs = 100*cm3Tom3;
    vs = system.branches.c0[1][end];
    if system.heart.av.zeta[n+1] == 0
        zs = 1e-5;
    else
        zs = system.heart.av.zeta[n+1];
    end

    # initial guess
    x0 = zeros(3);
    if system.branches.Q[1][n+1,1] == 0
        # x0[1] = 2*1e-6/system.branches.A[1][n+1,1]-system.branches.W2root;
        x0[1] = -system.branches.W2root + 1e-7;
    else
        x0[1] = system.branches.W1root;
        # x0[1] = (system.branches.W2root + 8*(sqrt(system.branches.beta[1][end]/
        #         (2*system.solverparams.rho))*system.branches.A[1][n+1,1]^0.25-
        #         system.branches.c0[1][end]));
    end
    if system.branches.Q[1][n+1,1] != 0
        x0[2] = (system.heart.lv.V[n+1] -
            system.solverparams.h*system.branches.Q[1][n+1,1]);
    elseif state == "opening" && system.branches.Q[1][n+1,1] == 0
        x0[2] = system.heart.lv.V[n+1]
        # x0[2] = (system.heart.lv.V[n+1] -
        #     system.solverparams.h*1e-6);
    end
    if system.heart.av.zeta[n+1] == 0
        x0[3] = zs;
    elseif state == "opening"
        x0[3] = system.heart.av.zeta[n+1] + 1e-5;
    elseif state == "closing" && abs(system.heart.lv.P[n+1]/mmHgToPa - system.branches.P[1][n+1,1]) < 20
        x0[3] = system.heart.av.zeta[n+1] - 2e-5;
    elseif state == "closing"
        x0[3] = system.heart.av.zeta[n+1];
    end

    # non-dimensionalize
    x0[1] = x0[1]/vs;
    x0[2] = x0[2]/Vs;
    x0[3] = x0[3]/zs;
    # println(x0)

    # assign function/Jacobian handles
    f = fav;
    J = Jav;

    # setup for iterations
    xx = x0;
    x = zeros(length(x0));
    N = 1;

    # max step size for line searches
    stpmax = 100*max(sqrt(norm(xx)),length(x0));
    # println(stpmax)

    # println(state)
    # Newton iterations
    while N <= system.solverparams.maxiter
        Plv = system.heart.lv.E[n+2]*(xx[2]*Vs-system.heart.lv.V0);
        Pao = system.branches.beta[1][end]*((2*system.solverparams.rho/
            system.branches.beta[1][end])*((xx[1]*vs-system.branches.W2root)/8+
            system.branches.c0[1][end])^2-sqrt(system.branches.A0[1][end]));
        if xx[3]*zs >= 0 && Plv > Pao
            state = "opening";
        else
            state = "closing";
            # println(state)
        end
        # println(xx)
        # determine Jacobian, check invertibility
        JJ = J(xx,system,n,state);
        D = diagm(maximum!(zeros(3),abs(JJ)).^-1);
        # println(JJ)
        # println(D)
        # println(D*JJ)
        # println(cond(D*JJ))
        # println(det(D*JJ))
        if abs(det(D*JJ)) < system.solverparams.epsJ
            println(JJ)
            println(D*JJ)
            print(det(D*JJ))
            error("Newton Jacobian is singular.");
        end
        # # update state vector
        # xn = xx - inv(D*JJ)*D*f(xx,system,n,state);
        # compute gradient of line search objective function
        fvec = D*f(xx,system,n,state);
        # println(f(xx,system,n,state))
        # println(fvec)
        fold = 0.5*dot(fvec,fvec);
        g = (D*JJ)'*fvec;
        # compute newton step
        # println(inv(D*JJ))
        s = -inv(D*JJ)*fvec;
        # line search to update state vector
        fn, xn, check = linesearch(xx,fold,g,s,stpmax,system,n,state);
        # println(xn[1]*vs)
        # check if sufficiently close to root
        JJ = J(xn,system,n,state);
        D = diagm(maximum!(zeros(3),abs(JJ)).^-1);
        fvec = D*f(xn,system,n,state);
        # println(fvec)
        if norm(fvec) <= system.solverparams.epsN*100
            x = xn;
            x[1] = x[1]*vs;
            x[2] = x[2]*Vs;
            x[3] = x[3]*zs;
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
    system.branches.A[1][n+2,1] = ((2*system.solverparams.rho
        /system.branches.beta[1][end])^2*((x[1]-system.branches.W2root)/8+
        system.branches.c0[1][end])^4);
    system.branches.Q[1][n+2,1] = 0.5*system.branches.A[1][n+2,1]*
        (x[1]+system.branches.W2root);
    system.heart.lv.V[n+2] = x[2];
    if x[3] > system.solverparams.epsN
        system.heart.av.zeta[n+2] = x[3];
    elseif x[3] < system.solverparams.epsN && state == "closing"
        system.heart.av.zeta[n+2] = 0;
    end

    system.branches.W1root = x[1];
end
