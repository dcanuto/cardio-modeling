function coupledistal!(system::CVSystem,n::Int64)
    j = Int8(1);
    omega = Float64(1);
    Qest = Float64;
    Vest = Float64;
    Pest = Float64;
    Aest = Float64;
    Uest = Float64;
    Qnew = Float64;
    W2 = Float64;
    iterr = Float64;

    # W1 at next time step
    endinvariants!(system,n);

    for i = 1:length(system.branches.ID)
        if isempty(system.branches.children[i])
            while j <= system.solverparams.maxiter
                if j == 1
                    Qest = system.branches.Q[i][n+1,system.solverparams.JL];
                end
                Vest = (system.branches.term[i].V[n+1,1]+
                    system.solverparams.h*(Qest-
                    system.branches.term[i].Q[n+1,1]));
                Pest = (1/system.branches.term[i].C[1]*
                    (Vest - system.branches.term[i].V0[1]));
                Aest = (Pest/system.branches.beta[i][end]+
                    sqrt(system.branches.A0[i][end]))^2;
                W2 = (system.branches.W1end[i] + 8*(
                    system.branches.c0[i][end]-
                    sqrt(system.branches.beta[i][end]/(2*
                    system.solverparams.rho))*Aest^0.25));
                Uest = 0.5*(system.branches.W1end[i] + W2);
                Qnew = (1-omega)*Qest + omega*Aest*Uest;
                iterr = abs(Qnew - Qest);
                if iterr <= system.solverparams.abserror
                    system.branches.A[i][n+2,system.solverparams.JL] = Aest;
                    system.branches.Q[i][n+2,system.solverparams.JL] = Qest;
                    system.branches.term[i].V[n+2,1] = (
                        system.branches.term[i].V[n+1,1]+
                        system.solverparams.h*(Qest-
                        system.branches.term[i].Q[n+1,1]));
                    system.branches.term[i].P[n+2,1] = Pest;
                    break
                else
                    Qest = Qnew;
                    j+=1;
                end
                if j == system.solverparams.maxiter
                    error("Distal coupling failed to converge.")
                end
            end
        end
    end
end
