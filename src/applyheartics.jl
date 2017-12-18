function applyheartics!(system::CVSystem,old=Dict("a"=>0),restart="no")
    if restart == "no"
        # volumes based on end-diastole
        system.heart.lv.V[1] = 120*cm3Tom3;
        # system.heart.rv.V[1] = 125*cm3Tom3;
        system.heart.rv.V[1] = 125*cm3Tom3;
        # system.heart.la.V[1] = 36*cm3Tom3;
        system.heart.la.V[1] = 36.5*cm3Tom3;
        system.heart.ra.V[1] = 32.7*cm3Tom3;
        # system.heart.ra.V[1] = 35.2*cm3Tom3;

        # aortic valve state
        system.heart.av.zeta[1] = 0;
    elseif restart == "yes"
        lv = old["lv"];
        rv = old["rv"];
        la = old["la"];
        ra = old["ra"];
        av = old["av"];
        system.heart.lv.V[1] = lv["V"][end];
        system.heart.rv.V[1] = rv["V"][end];
        system.heart.la.V[1] = la["V"][end];
        system.heart.ra.V[1] = ra["V"][end];
        system.heart.av.zeta[1] = av["zeta"][end];
    end

    # associated pressures
    elastancemodel!(system,0);
end
