function updatevc!(system::CVSystem,n::Int64)
    ivcflow = Float64(0);
    svcflow = Float64(0);
    haveupper = Int8(0);

    for i = 1:length(system.branches.ID)
        if isempty(system.branches.children[i])
            if system.branches.group[i] == "lower"
                ivcflow += system.branches.term[i].Q[n+1,5];
            else
                svcflow += system.branches.term[i].Q[n+1,5];
                haveupper = 1;
            end
        end
    end

    if haveupper == 1
        system.ivc.V[n+2] = (system.ivc.V[n+1] + system.solverparams.h*(
            ivcflow - system.ivc.Q[n+1]));
        system.ivc.Q[n+2] = ((1 - system.ivc.R/system.ivc.L*
            system.solverparams.h)*system.ivc.Q[n+1]+
            system.solverparams.h/system.ivc.L*(system.ivc.P[n+1]-
            system.heart.ra.P[n+1]));
        system.ivc.P[n+2] = (1/system.ivc.C*(system.ivc.V[n+2]-
            system.ivc.V0));
        system.svc.V[n+2] = (system.svc.V[n+1] + system.solverparams.h*(
            svcflow - system.svc.Q[n+1]));
        system.svc.Q[n+2] = ((1 - system.svc.R/system.svc.L*
            system.solverparams.h)*system.svc.Q[n+1]+
            system.solverparams.h/system.svc.L*(system.svc.P[n+1]-
            system.heart.ra.P[n+1]));
        system.svc.P[n+2] = (1/system.svc.C*(system.svc.V[n+2]-
            system.svc.V0));
    else
        system.ivc.V[n+2] = (system.ivc.V[n+1] + system.solverparams.h*(
            ivcflow - system.ivc.Q[n+1]));
        system.ivc.Q[n+2] = ((1 - system.ivc.R/system.ivc.L*
            system.solverparams.h)*system.ivc.Q[n+1]+
            system.solverparams.h/system.ivc.L*(system.ivc.P[n+1]-
            system.heart.ra.P[n+1]));
        system.ivc.P[n+2] = (1/system.ivc.C*(system.ivc.V[n+2]-
            system.ivc.V0));
    end
end
