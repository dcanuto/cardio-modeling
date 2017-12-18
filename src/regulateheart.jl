function regulateheart!(system::CVSystem,n::Int64)
    # change activation function parameters based on new HR
    push!(system.heart.activation.th,1/system.cns.H[n+2]*minTos);
    system.heart.activation.tau1 = 0.269*system.heart.activation.th[end];
    system.heart.activation.tau2 = 0.452*system.heart.activation.th[end];
    t = linspace(0,system.heart.activation.th[end],10000);
    g1 = (t/system.heart.activation.tau1).^system.heart.activation.m1;
    g2 = (t/system.heart.activation.tau2).^system.heart.activation.m2;
    push!(system.heart.activation.k,maximum((g1./(1+g1)).*(1./(1+g2)))^-1)

    # change max. ventriclular elastance
    push!(system.heart.lv.Emax,system.cns.Emaxlv[n+2]);
    push!(system.heart.rv.Emax,system.cns.Emaxrv[n+2]);
end
