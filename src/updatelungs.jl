function updatelungs!(system::CVSystem,n::Int64)
    # 1st arterial compartment
    system.lungs.Va[n+2,1] = (system.lungs.Va[n+1,1] + system.solverparams.h*
        (system.heart.rv.Q[n+1] - system.lungs.Qa[n+1,1]));
    system.lungs.Qa[n+2,1] = ((1 - system.lungs.Ra[1]/system.lungs.La[1]*
        system.solverparams.h)*system.lungs.Qa[n+1,1]+
        system.solverparams.h/system.lungs.La[1]*(system.lungs.Pa[n+1,1]-
        system.lungs.Pa[n+1,2]));
    system.lungs.Pa[n+2,1] = (1/system.lungs.Ca[1]*(system.lungs.Va[n+2,1]-
        system.lungs.V0[1]));

    # 2nd arterial compartment V, P
    system.lungs.Va[n+2,2] = (system.lungs.Va[n+1,2] + system.solverparams.h*
        (system.lungs.Qa[n+1,1] - system.lungs.Qa[n+1,2]));
    system.lungs.Pa[n+2,2] = (1/system.lungs.Ca[2]*(system.lungs.Va[n+2,2]-
        system.lungs.V0[2]));

    # 3rd arterial compartment V, P
    system.lungs.Va[n+2,3] = (system.lungs.Va[n+1,3] + system.solverparams.h*
        (system.lungs.Qa[n+1,2] - system.lungs.Qa[n+1,3]));
    system.lungs.Pa[n+2,3] = (1/system.lungs.Ca[3]*(system.lungs.Va[n+2,3]-
        system.lungs.V0[3]));

    # 1st venous compartment V, P
    system.lungs.Vv[n+2,1] = (system.lungs.Vv[n+1,1] + system.solverparams.h*
        (system.lungs.Qa[n+1,3] - system.lungs.Qv[n+1,1]));
    system.lungs.Pv[n+2,1] = (1/system.lungs.Cv[1]*(system.lungs.Vv[n+2,1]-
        system.lungs.V0[4]));

    # 2nd venous compartment
    system.lungs.Vv[n+2,2] = (system.lungs.Vv[n+1,2] + system.solverparams.h*
        (system.lungs.Qv[n+1,1] - system.lungs.Qv[n+1,2]));
    system.lungs.Qv[n+2,2] = ((1 - system.lungs.Rv[2]/system.lungs.Lv[1]*
        system.solverparams.h)*system.lungs.Qv[n+1,2]+
        system.solverparams.h/system.lungs.Lv[1]*(system.lungs.Pv[n+1,2]-
        system.heart.la.P[n+1]));
    system.lungs.Pv[n+2,2] = (1/system.lungs.Cv[2]*(system.lungs.Vv[n+2,2]-
        system.lungs.V0[5]));

    # remaining Q's
    system.lungs.Qa[n+2,2] = ((system.lungs.Pa[n+2,2] - system.lungs.Pa[n+2,3])/
        system.lungs.Ra[2]);
    system.lungs.Qa[n+2,3] = ((system.lungs.Pa[n+2,3] - system.lungs.Pv[n+2,1])/
        system.lungs.Ra[3]);
    system.lungs.Qv[n+2,1] = ((system.lungs.Pv[n+2,1] - system.lungs.Pv[n+2,2])/
        system.lungs.Rv[1]);

    # root pulmonary artery P
    system.lungs.Pp[n+2] = (system.lungs.Rp*system.heart.rv.Q[n+2] +
        system.lungs.Pa[n+2,1]);
end
