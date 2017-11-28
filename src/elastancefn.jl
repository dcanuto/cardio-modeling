function elastancefn!(system::CVSystem,n::Int64)
    # time within heart cycle
    if system.solverparams.numbeats > 0
        tp = system.t[n+1] - sum(system.heart.activation.th[1:end-1]);
    else
        tp = system.t[n+1];
    end

    # heart activation function
    if tp < system.heart.activation.tce[system.solverparams.numbeats+1]
        phi = (system.heart.activation.a*sin(π*tp/
            system.heart.activation.tce[system.solverparams.numbeats+1])-
            system.heart.activation.b*sin(2*π*tp/
            system.heart.activation.tce[system.solverparams.numbeats+1]));
    else
        phi = 0;
    end

    # ventricular elastance
    system.heart.lv.E[n+1] = (system.heart.lv.Emin*(1-phi) +
        system.heart.lv.Emax[system.solverparams.numbeats+1]*phi);
    system.heart.rv.E[n+1] = (system.heart.rv.Emin*(1-phi) +
        system.heart.rv.Emax[system.solverparams.numbeats+1]*phi);
end
