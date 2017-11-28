function elastancemodel!(system::CVSystem,n::Int64)
    # ventricular elastance
    elastancefn!(system,n)

    # update chamber pressures
    system.heart.lv.P[n+1] = system.heart.lv.E[n+1]*(
        system.heart.lv.V[n+1] - system.heart.lv.V0);
    system.heart.rv.P[n+1] = system.heart.rv.E[n+1]*(
        system.heart.rv.V[n+1] - system.heart.rv.V0);

    system.heart.la.P[n+1] = system.heart.la.E*(
        system.heart.la.V[n+1] - system.heart.la.V0);
    system.heart.ra.P[n+1] = system.heart.ra.E*(
        system.heart.ra.V[n+1] - system.heart.ra.V0);
end
