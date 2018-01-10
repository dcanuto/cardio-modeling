function linedist(xold::Vector{Float64},fold::Float64,
    g::Vector{Float64},p::Vector{Float64},stpmax::Float64,
    system::CVSystem,n::Int64,ID::Int64)
    alpha = 1e-4;
    tolx = 1e-7;
    check = false;

    f = fdist;
    J = Jdist;

    # scale if attempted step too large
    np = norm(p);
    # println(p)
    if np > stpmax
        p = p*stpmax/np;
    end
    # println(p)

    slope = dot(g,p);
    if slope >= 0
        print(g)
        print(p)
        error("Roundoff problem in line search.")
    end

    # compute minimum step fraction
    test = 0.;
    for i = 1:length(xold)
        temp = abs(p[i])/max(abs(xold[i]),1.);
        if temp > test
            test = temp;
        end
    end
    alamin = tolx/test;

    omega = 1;
    alam = omega; # start w/ desired fraction of Newton step
    N = 1;

    while N < system.solverparams.maxiter
        x = xold+alam*p;
        JJ = J(x,system,n,ID);
        # println(JJ)
        D = diagm(maximum!(zeros(length(x)),abs(JJ)).^-1);
        fvec = D*f(x,system,n,ID);
        # println(D)
        # println(fvec)
        fn = 0.5*dot(fvec,fvec);
        if alam < alamin
            x = xold;
            # println(alam)
            # println(alamin)
            # println(fn)
            # println(fold+alpha*alam*slope)
            # println("Δx converged in distal line search. Verify in Newton loop.")
            check = true;
            return fn,x,check
        elseif fn < fold+alpha*alam*slope
            # println(fn)
            # println(fold)
            # println(alpha*alam*slope)
            return fn,x,check
        else
            if alam == omega
                tmplam = -slope/(2*(fn-fold-slope));
            else
                rhs1 = fn-fold-alam*slope
                rhs2 = fn2-fold-alam2*slope
                a = (rhs1/alam^2-rhs2/alam2^2)/(alam-alam2);
                b = (-alam2*rhs1/alam^2+alam*rhs2/alam2^2)/(alam-alam2);
                if a == 0
                    tmplam = -slope/(2*b);
                else
                    disc = b^2-3*a*slope;
                    if disc < 0
                        tmplam = 0.5*alam;
                    elseif b <= 0
                        tmplam = (-b+sqrt(disc)/(3*a));
                    else
                        tmplam = -slope/(b+sqrt(disc));
                    end
                end
                if tmplam > 0.5*alam
                    tmplam = 0.5*alam;
                end
            end
        end
        alam2 = alam;
        fn2 = fn;
        alam = max(tmplam,0.1*alam);
        N+=1;
        if N == system.solverparams.maxiter
            println(xn)
            println(f(xn,system,n,ID))
            error("Distal line search iteration failed to converge.");
        end
    end
end
