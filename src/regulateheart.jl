function regulateheart!(system::CVSystem,n::Int64)
    # change HR, time for const. elastance
    push!(system.heart.activation.th,1/system.cns.H[n+2]*minTos);
    push!(system.heart.activation.tce,system.heart.activation.k0+
        system.heart.activation.k1*system.heart.activation.th[end]);

    # change max. ventriclular elastance
    push!(system.heart.lv.Emax,system.cns.Emaxlv[n+2]);
    push!(system.heart.rv.Emax,system.cns.Emaxrv[n+2]);
end
