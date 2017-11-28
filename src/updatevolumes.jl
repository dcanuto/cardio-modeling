function updatevolumes!(system::CVSystem,n::Int64)
    system.arterialvolume = 0;
    system.peripheralvolume = 0;
    system.vcvolume = 0;
    system.heartvolume = 0;
    system.lungvolume = 0;

    # 1D arterial domain and peripheral compartments
    for i = 1:length(system.branches.ID)
        x = convert(Array{Float64,1},
            linspace(0,system.branches.lengthincm[i]*cmTom,
            system.solverparams.JL));
        system.arterialvolume += integrate(x,system.branches.A[i][n+1,:]);
        if isempty(system.branches.children[i])
            system.arterialvolume += system.branches.term[i].V[n+1,1];
            system.peripheralvolume += sum(system.branches.term[i].V[n+1,2:5]);
        end
    end
     # vena cava, heart, and lung compartments
    system.vcvolume += (system.svc.V[n+1] + system.ivc.V[n+1]);
    system.heartvolume += (system.heart.rv.V[n+1] + system.heart.ra.V[n+1] +
        system.heart.lv.V[n+1] + system.heart.la.V[n+1]);
    system.lungvolume += (sum(system.lungs.Va[n+1,:]) + sum(system.lungs.Vv[n+1,:]));

    if n == 0
        system.initialvolume = (system.arterialvolume +
            system.peripheralvolume + system.vcvolume + system.heartvolume +
            system.lungvolume);
    else
        system.finalvolume = (system.arterialvolume +
            system.peripheralvolume + system.vcvolume + system.heartvolume +
            system.lungvolume);
    end
end
