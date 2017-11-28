function updatela!(system::CVSystem,n::Int64)
    # mitral valve state
    if system.heart.la.P[n+1] > system.heart.lv.P[n+1]
        MV = "open"
    else
        MV = "closed"
    end

    # left atrial V, Q (P updated elsewhere)
    system.heart.la.V[n+2] = (system.heart.la.V[n+1]+
        system.solverparams.h*(system.lungs.Qv[n+1,2]-
        system.heart.la.Q[n+1]));
    if MV == "open"
        system.heart.la.Q[n+2] = ((1 - system.heart.la.R/
            system.heart.la.L*system.solverparams.h)*
            system.heart.la.Q[n+1] + system.solverparams.h/
            system.heart.la.L*(system.heart.la.P[n+1]-
            system.heart.lv.P[n+1]));
    else
        system.heart.la.Q[n+2] = 0;
    end
end
