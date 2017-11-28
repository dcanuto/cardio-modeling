function applyheartics!(system::CVSystem)
    # volumes based on end-diastole
    system.heart.lv.V[1] = 120*cm3Tom3;
    system.heart.rv.V[1] = 125*cm3Tom3;
    system.heart.la.V[1] = 36*cm3Tom3;
    system.heart.ra.V[1] = 32.7*cm3Tom3;

    # associated pressures
    elastancemodel!(system,0);
end
