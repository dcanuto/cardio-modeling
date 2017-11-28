function updaterh!(system::CVSystem,n::Int64)
    # pulmonary, tricuspid valve states
    if system.heart.ra.P[n+1] > system.heart.rv.P[n+1]
        TV = "open";
    else
        TV = "closed";
    end

    if (system.heart.rv.P[n+1] > system.lungs.Pp[n+1] ||
        system.heart.rv.Q[n+1] > 0)
        PV = "open";
    else
        PV = "closed";
    end

    # right atrial V, Q (P updated elsewhere)
    system.heart.ra.V[n+2] = (system.heart.ra.V[n+1]+
        system.solverparams.h*(system.ivc.Q[n+1] + system.svc.Q[n+1]-
        system.heart.ra.Q[n+1]));
    if TV == "open"
        system.heart.ra.Q[n+2] = ((1 - system.heart.ra.R/system.heart.ra.L*
            system.solverparams.h)*system.heart.ra.Q[n+1]+
            system.solverparams.h/system.heart.ra.L*(system.heart.ra.P[n+1]-
            system.heart.rv.P[n+1]));
    else
        system.heart.ra.Q[n+2] = 0;
    end

    # right ventricular V, Q
    system.heart.rv.V[n+2] = (system.heart.rv.V[n+1]+
        system.solverparams.h*(system.heart.ra.Q[n+1]-
        system.heart.rv.Q[n+1]));
    if PV == "open"
        system.heart.rv.Q[n+2] = (system.heart.rv.Q[n+1]+
            system.solverparams.h/system.heart.rv.L*(system.heart.rv.P[n+1]-
            system.lungs.Pp[n+1]));
    else
        system.heart.rv.Q[n+2] = 0;
    end
end
