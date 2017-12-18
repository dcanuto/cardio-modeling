function updatediscretization!(system::CVSystem)
    # determine time shift
    system.solverparams.tshift = system.t[end]-
        sum(system.heart.activation.th[1:end-1]);

    # update total number of time steps
    ntoadd = ceil((system.heart.activation.th[end]-
        system.solverparams.tshift)/system.solverparams.h);
    system.solverparams.numsteps+=ntoadd;

    # update discrete times
    ttoadd = [0:1:ntoadd;]*system.solverparams.h + system.t[end];
    append!(system.t,ttoadd[2:end]);

    # allocate additional space for solution variables
    for i = 1:length(system.branches.ID)
        system.branches.A[i] = [system.branches.A[i];
            zeros(length(ttoadd[2:end]),system.solverparams.JL)];
        system.branches.Q[i] = [system.branches.Q[i];
            zeros(length(ttoadd[2:end]),system.solverparams.JL)];
        system.branches.P[i] = [system.branches.P[i];
            zeros(length(ttoadd[2:end]),system.solverparams.JL)];
        if isempty(system.branches.children[i])
            system.branches.term[i].P = [system.branches.term[i].P;
                zeros(length(ttoadd[2:end]),5)];
            system.branches.term[i].V = [system.branches.term[i].V;
                zeros(length(ttoadd[2:end]),5)];
            system.branches.term[i].Q = [system.branches.term[i].Q;
                zeros(length(ttoadd[2:end]),5)];
        end
    end

    append!(system.svc.P,zeros(length(ttoadd[2:end])));
    append!(system.svc.V,zeros(length(ttoadd[2:end])));
    append!(system.svc.Q,zeros(length(ttoadd[2:end])));
    append!(system.ivc.P,zeros(length(ttoadd[2:end])));
    append!(system.ivc.V,zeros(length(ttoadd[2:end])));
    append!(system.ivc.Q,zeros(length(ttoadd[2:end])));

    append!(system.heart.lv.V,zeros(length(ttoadd[2:end])));
    append!(system.heart.lv.P,zeros(length(ttoadd[2:end])));
    append!(system.heart.lv.E,zeros(length(ttoadd[2:end])));

    append!(system.heart.av.zeta,zeros(length(ttoadd[2:end])));

    append!(system.heart.rv.V,zeros(length(ttoadd[2:end])));
    append!(system.heart.rv.P,zeros(length(ttoadd[2:end])));
    append!(system.heart.rv.E,zeros(length(ttoadd[2:end])));
    append!(system.heart.rv.Q,zeros(length(ttoadd[2:end])));

    append!(system.heart.la.V,zeros(length(ttoadd[2:end])));
    append!(system.heart.la.P,zeros(length(ttoadd[2:end])));
    append!(system.heart.la.Q,zeros(length(ttoadd[2:end])));

    append!(system.heart.ra.V,zeros(length(ttoadd[2:end])));
    append!(system.heart.ra.P,zeros(length(ttoadd[2:end])));
    append!(system.heart.ra.Q,zeros(length(ttoadd[2:end])));

    append!(system.lungs.Pp,zeros(length(ttoadd[2:end])));
    system.lungs.Pa = [system.lungs.Pa;zeros(length(ttoadd[2:end]),3)];
    system.lungs.Va = [system.lungs.Va;zeros(length(ttoadd[2:end]),3)];
    system.lungs.Qa = [system.lungs.Qa;zeros(length(ttoadd[2:end]),3)];
    system.lungs.Pv = [system.lungs.Pv;zeros(length(ttoadd[2:end]),2)];
    system.lungs.Vv = [system.lungs.Vv;zeros(length(ttoadd[2:end]),2)];
    system.lungs.Qv = [system.lungs.Qv;zeros(length(ttoadd[2:end]),2)];

    append!(system.cns.H,zeros(length(ttoadd[2:end])));
    append!(system.cns.Emaxlv,zeros(length(ttoadd[2:end])));
    append!(system.cns.Emaxrv,zeros(length(ttoadd[2:end])));
    append!(system.cns.R2L,zeros(length(ttoadd[2:end])));
    append!(system.cns.R3L,zeros(length(ttoadd[2:end])));
    append!(system.cns.R2U,zeros(length(ttoadd[2:end])));
    append!(system.cns.R3U,zeros(length(ttoadd[2:end])));
    append!(system.cns.C4L,zeros(length(ttoadd[2:end])));
    append!(system.cns.C5L,zeros(length(ttoadd[2:end])));
    append!(system.cns.C4U,zeros(length(ttoadd[2:end])));
    append!(system.cns.C5U,zeros(length(ttoadd[2:end])));
    append!(system.cns.V4L,zeros(length(ttoadd[2:end])));
    append!(system.cns.V5L,zeros(length(ttoadd[2:end])));
    append!(system.cns.V4U,zeros(length(ttoadd[2:end])));
    append!(system.cns.V5U,zeros(length(ttoadd[2:end])));

end
