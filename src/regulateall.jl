function regulateall!(system::CVSystem,n::Int64)
    # changes to HR, ventricular contractility
    system.cns.H[n+2] = ((1 - system.solverparams.h/system.cns.tauH)*
        system.cns.H[n+1]/minTos + system.solverparams.h/system.cns.tauH*
        (system.cns.alphaH*system.cns.ns[system.solverparams.numbeats+1]-
        system.cns.betaH*system.cns.np[system.solverparams.numbeats+1]+
        system.cns.gammaH))*minTos;
    system.cns.Emaxlv[n+2] = (1 - system.solverparams.h/system.cns.tauEmax)*
        system.cns.Emaxlv[n+1] + system.solverparams.h/system.cns.tauEmax*(
        system.cns.alphalv*system.cns.ns[system.solverparams.numbeats+1]+
        system.cns.gammalv);
    system.cns.Emaxrv[n+2] = (1 - system.solverparams.h/system.cns.tauEmax)*
        system.cns.Emaxrv[n+1] + system.solverparams.h/system.cns.tauEmax*(
        system.cns.alpharv*system.cns.ns[system.solverparams.numbeats+1]+
        system.cns.gammarv);

    # changes to peripheral impedance
    system.cns.R2L[n+2] = (1 - system.solverparams.h/system.cns.tauR)*
        system.cns.R2L[n+1] + system.solverparams.h/system.cns.tauR*(
        system.cns.alphaR2L*system.cns.ns[system.solverparams.numbeats+1]+
        system.cns.gammaR2L);
    system.cns.R3L[n+2] = (1 - system.solverparams.h/system.cns.tauR)*
        system.cns.R3L[n+1] + system.solverparams.h/system.cns.tauR*(
        system.cns.alphaR3L*system.cns.ns[system.solverparams.numbeats+1]+
        system.cns.gammaR3L);
    system.cns.R2U[n+2] = (1 - system.solverparams.h/system.cns.tauR)*
        system.cns.R2U[n+1] + system.solverparams.h/system.cns.tauR*(
        system.cns.alphaR2U*system.cns.ns[system.solverparams.numbeats+1]+
        system.cns.gammaR2U);
    system.cns.R3U[n+2] = (1 - system.solverparams.h/system.cns.tauR)*
        system.cns.R3U[n+1] + system.solverparams.h/system.cns.tauR*(
        system.cns.alphaR3U*system.cns.ns[system.solverparams.numbeats+1]+
        system.cns.gammaR3U);
    system.cns.C4L[n+2] = (1 - system.solverparams.h/system.cns.tauC)*
        system.cns.C4L[n+1] + system.solverparams.h/system.cns.tauC*(
        -system.cns.alphaC4L*system.cns.ns[system.solverparams.numbeats+1]+
        system.cns.gammaC4L);
    system.cns.C5L[n+2] = (1 - system.solverparams.h/system.cns.tauC)*
        system.cns.C5L[n+1] + system.solverparams.h/system.cns.tauC*(
        -system.cns.alphaC5L*system.cns.ns[system.solverparams.numbeats+1]+
        system.cns.gammaC5L);
    system.cns.C4U[n+2] = (1 - system.solverparams.h/system.cns.tauC)*
        system.cns.C4U[n+1] + system.solverparams.h/system.cns.tauC*(
        -system.cns.alphaC4U*system.cns.ns[system.solverparams.numbeats+1]+
        system.cns.gammaC4U);
    system.cns.C5U[n+2] = (1 - system.solverparams.h/system.cns.tauC)*
        system.cns.C5U[n+1] + system.solverparams.h/system.cns.tauC*(
        -system.cns.alphaC5U*system.cns.ns[system.solverparams.numbeats+1]+
        system.cns.gammaC5U);
    system.cns.V4L[n+2] = (1 - system.solverparams.h/system.cns.tauC)*
        system.cns.V4L[n+1] + system.solverparams.h/system.cns.tauC*(
        -system.cns.alphaV4L*system.cns.ns[system.solverparams.numbeats+1]+
        system.cns.gammaV4L);
    system.cns.V5L[n+2] = (1 - system.solverparams.h/system.cns.tauC)*
        system.cns.V5L[n+1] + system.solverparams.h/system.cns.tauC*(
        -system.cns.alphaV5L*system.cns.ns[system.solverparams.numbeats+1]+
        system.cns.gammaV5L);
    system.cns.V4U[n+2] = (1 - system.solverparams.h/system.cns.tauC)*
        system.cns.V4U[n+1] + system.solverparams.h/system.cns.tauC*(
        -system.cns.alphaV4U*system.cns.ns[system.solverparams.numbeats+1]+
        system.cns.gammaV4U);
    system.cns.V5U[n+2] = (1 - system.solverparams.h/system.cns.tauC)*
        system.cns.V5U[n+1] + system.solverparams.h/system.cns.tauC*(
        -system.cns.alphaV5U*system.cns.ns[system.solverparams.numbeats+1]+
        system.cns.gammaV5U);

    regulateperiphery!(system,n);

    # check for start of new cardiac cycle
    oldnumbeats = system.solverparams.numbeats;
    setnumbeats!(system,n);
    if oldnumbeats < system.solverparams.numbeats
        # baroreflex pressure
        reflexpressure!(system,n);
        # SNS/PSNS activations
        cnsactivations!(system,n);
        regulateheart!(system,n);
        # update time discretization
        if system.solverparams.numbeats < system.solverparams.numbeatstotal
            updatediscretization!(system);
        end
    end
end
