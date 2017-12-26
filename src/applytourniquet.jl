function applytourniquet!(system::CVSystem,n::Int64)
    if system.hemo.Vloss >= system.hemo.totalloss
        for i = 1:length(system.branches.ID)
            tb = system.branches.beta[i][end];
            ta = system.branches.A0[i][end];
            tc = system.branches.c0[i][end]
            for j = 1:length(system.hemo.tID)
                if i == system.hemo.tID[j]
                    if system.branches.beta[i][end] < system.hemo.bmax[j]
                        tb = (system.branches.beta[i][end]+
                            (system.hemo.bmax[j]-system.hemo.bmin[j])/system.hemo.ttn*
                            system.solverparams.h);
                    end
                    if system.branches.A0[i][end] > system.hemo.Amin[j]
                        ta = (system.branches.A0[i][end]+
                            (system.hemo.Amin[j]-system.hemo.Amax[j])/system.hemo.ttn*
                            system.solverparams.h);
                        tc = sqrt(system.branches.beta[i][end]/
                            (2*system.solverparams.rho))*system.branches.A0[i][end]^0.25;
                    end
                end
            end
            push!(system.branches.beta[i],tb)
            push!(system.branches.A0[i],ta);
            push!(system.branches.c0[i],tc);
        end
    end
end
