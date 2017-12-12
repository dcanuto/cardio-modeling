function Jav(x::Vector{Float64},system::CVSystem,n::Int64,state::String)
    # non-dimensionalizing parameters
    Vs = 100*cm3Tom3;
    vs = system.branches.c0[1][end];
    # vs = 1;
    ts = system.heart.av.leff/vs;
    if system.heart.av.zeta[n+1] == 0
        zs = 1e-5;
    else
        zs = system.heart.av.zeta[n+1];
    end

    J11 = (-2*ts*vs^5/Vs*(system.solverparams.rho/system.branches.beta[1][end])^2*
        (((x[1]-system.branches.W2root/vs)/8+system.branches.c0[1][end]/vs)^3*
        ((x[1]+system.branches.W2root/vs)/2)+((x[1]-system.branches.W2root/vs)/8+
        system.branches.c0[1][end]/vs)^4));
    J12 = -1*ts/system.solverparams.h;
    J13 = 0;
    J1 = [J11 J12 J13];

    if state == "opening"
        J21 = (0.5*system.solverparams.h*(zs^-1-x[3])*system.solverparams.rho*vs^2*
            system.heart.av.Kvo*((x[1]-system.branches.W2root/vs)/8+
            system.branches.c0[1][end]/vs));
        J22 = Vs*(-system.solverparams.h*(zs^-1-x[3])*system.heart.av.Kvo*
            system.heart.lv.E[n+2]);
        J23 = (1+system.solverparams.h*system.heart.av.Kvo*(system.heart.lv.E[n+2]*
            Vs*(x[2]-system.heart.lv.V0/Vs)-system.branches.beta[1][end]*
            ((2*system.solverparams.rho/system.branches.beta[1][end])*vs^2*
            ((x[1]-system.branches.W2root/vs)/8+system.branches.c0[1][end]/vs)^2-
            sqrt(system.branches.A0[1][end]))));
    else
        J21 = (0.5*system.solverparams.h*(x[3])*system.solverparams.rho*vs^2*
            system.heart.av.Kvc*((x[1]-system.branches.W2root/vs)/8+
            system.branches.c0[1][end]/vs));
        J22 = Vs*(-system.solverparams.h*(x[3])*system.heart.av.Kvc*
            system.heart.lv.E[n+2]);
        J23 = (1+system.solverparams.h*system.heart.av.Kvc*(system.heart.lv.E[n+2]*
            Vs*(x[2]-system.heart.lv.V0/Vs)-system.branches.beta[1][end]*
            ((2*system.solverparams.rho/system.branches.beta[1][end])*vs^2*
            ((x[1]-system.branches.W2root/vs)/8+system.branches.c0[1][end]/vs)^2-
            sqrt(system.branches.A0[1][end]))));
    end
    J2 = [J21 J22 J23];

    J31 = 0;
    # J31 = ts*system.solverparams.h*system.heart.av.Aann*vs^2/(2*
    #     system.heart.av.leff*Vs*zs)*x[3]^2*((x[1]-system.branches.W2root/vs)/8+
    #     system.branches.c0[1][end]);
    # if x[2] > system.heart.lv.V[n+1]
    #     J32 = (-x[3]^3*ts/system.solverparams.h+Vs*ts/(2*system.heart.av.Aann*
    #         system.heart.av.leff*system.solverparams.h*zs^3)*
    #         ((system.heart.lv.V[n+1]/Vs-x[2])-abs(system.heart.lv.V[n+1]/Vs-x[2]))+
    #         x[3]^2*system.heart.av.Ks*system.heart.lv.E[n+2]*Vs*ts*
    #         system.heart.av.Aann/(system.solverparams.rho*system.heart.av.leff*zs)*
    #         (system.heart.lv.V[n+1]/Vs+system.heart.lv.V0/Vs-2*x[2]));
    #     # println(-x[3]^3*ts/system.solverparams.h)
    #     # println(Vs*ts/(2*system.heart.av.Aann*
    #     #     system.heart.av.leff*system.solverparams.h*zs^3)*
    #     #     ((system.heart.lv.V[n+1]/Vs-x[2])-abs(system.heart.lv.V[n+1]/Vs-x[2])))
    #     # println(x[3]^2*system.heart.av.Ks*system.heart.lv.E[n+2]*Vs*ts*
    #     #     system.heart.av.Aann/(system.solverparams.rho*system.heart.av.leff*zs)*
    #     #     (system.heart.lv.V[n+1]/Vs+system.heart.lv.V0/Vs-2*x[2]))
    #     # J32 = (-x[3]^3/system.solverparams.h*ts-ts*system.solverparams.h*
    #     #     system.heart.av.Aann*system.heart.lv.E[n+2]/(system.solverparams.rho*
    #     #     system.heart.av.leff*Vs*zs)*x[3]^2+ts*Vs/(2*system.heart.av.Aann*
    #     #     system.heart.av.leff*system.solverparams.h*zs^3)*((system.heart.lv.V[n+1]/Vs-
    #     #     x[2])-abs(system.heart.lv.V[n+1]/Vs-x[2]))+ts*Vs*system.heart.av.Aann*
    #     #     system.heart.av.Ks*system.heart.lv.E[n+2]/(system.solverparams.rho*
    #     #     system.heart.av.leff*zs)*x[3]^2*(system.heart.lv.V[n+1]/Vs+
    #     #     system.heart.lv.V0/Vs-2*x[2]));
    # else
    #     J32 = (-x[3]^3*ts/system.solverparams.h+Vs*ts/(2*system.heart.av.Aann*
    #         system.heart.av.leff*system.solverparams.h*zs^3)*
    #         ((-system.heart.lv.V[n+1]/Vs+x[2])-abs(system.heart.lv.V[n+1]/Vs-x[2]))+
    #         x[3]^2*system.heart.av.Ks*system.heart.lv.E[n+2]*Vs*ts*
    #         system.heart.av.Aann/(system.solverparams.rho*system.heart.av.leff*zs)*
    #         (system.heart.lv.V[n+1]/Vs+system.heart.lv.V0/Vs-2*x[2]));
    #     # J32 = (-x[3]^3/system.solverparams.h*ts-ts*system.solverparams.h*
    #     #     system.heart.av.Aann*system.heart.lv.E[n+2]/(system.solverparams.rho*
    #     #     system.heart.av.leff*Vs*zs)*x[3]^2+ts*Vs/(2*system.heart.av.Aann*
    #     #     system.heart.av.leff*system.solverparams.h*zs^3)*((-system.heart.lv.V[n+1]/Vs+
    #     #     x[2])-abs(system.heart.lv.V[n+1]/Vs-x[2]))+ts*Vs*system.heart.av.Aajnn*
    #     #     system.heart.av.Ks*system.heart.lv.E[n+2]/(system.solverparams.rho*
    #     #     system.heart.av.leff*zs)*x[3]^2*(system.heart.lv.V[n+1]/Vs+
    #     #     system.heart.lv.V0/Vs-2*x[2]));
    # end
    if state == "closing" && x[2] > system.heart.lv.V[n+1]
        J32 = (-x[3]^3*ts/system.solverparams.h+Vs*ts/(2*system.heart.av.Aann*
            system.heart.av.leff*system.solverparams.h*zs^3)*
            ((system.heart.lv.V[n+1]/Vs-x[2])-abs(system.heart.lv.V[n+1]/Vs-x[2]))+
            x[3]^2*system.heart.av.Ks*system.heart.lv.E[n+2]*Vs*ts*
            system.heart.av.Aann/(system.solverparams.rho*system.heart.av.leff*zs)*
            (system.heart.lv.V[n+1]/Vs+system.heart.lv.V0/Vs-2*x[2]));
        # println(-x[3]^3*ts/system.solverparams.h)
        # println(Vs*ts/(2*system.heart.av.Aann*
        #     system.heart.av.leff*system.solverparams.h*zs^3)*
        #     ((system.heart.lv.V[n+1]/Vs-x[2])-abs(system.heart.lv.V[n+1]/Vs-x[2])))
        # println(x[3]^2*system.heart.av.Ks*system.heart.lv.E[n+2]*Vs*ts*
        #     system.heart.av.Aann/(system.solverparams.rho*system.heart.av.leff*zs)*
        #     (system.heart.lv.V[n+1]/Vs+system.heart.lv.V0/Vs-2*x[2]))
        # J32 = (-x[3]^3/system.solverparams.h*ts-ts*system.solverparams.h*
        #     system.heart.av.Aann*system.heart.lv.E[n+2]/(system.solverparams.rho*
        #     system.heart.av.leff*Vs*zs)*x[3]^2+ts*Vs/(2*system.heart.av.Aann*
        #     system.heart.av.leff*system.solverparams.h*zs^3)*((system.heart.lv.V[n+1]/Vs-
        #     x[2])-abs(system.heart.lv.V[n+1]/Vs-x[2]))+ts*Vs*system.heart.av.Aann*
        #     system.heart.av.Ks*system.heart.lv.E[n+2]/(system.solverparams.rho*
        #     system.heart.av.leff*zs)*x[3]^2*(system.heart.lv.V[n+1]/Vs+
        #     system.heart.lv.V0/Vs-2*x[2]));
    elseif state == "opening"
        J32 = (-x[3]^3*ts/system.solverparams.h+Vs*ts/(2*system.heart.av.Aann*
            system.heart.av.leff*system.solverparams.h*zs^3)*
            ((-system.heart.lv.V[n+1]/Vs+x[2])-abs(system.heart.lv.V[n+1]/Vs-x[2]))+
            x[3]^2*system.heart.av.Ks*system.heart.lv.E[n+2]*Vs*ts*
            system.heart.av.Aann/(system.solverparams.rho*system.heart.av.leff*zs)*
            (system.heart.lv.V[n+1]/Vs+system.heart.lv.V0/Vs-2*x[2]));
        # J32 = (-x[3]^3/system.solverparams.h*ts-ts*system.solverparams.h*
        #     system.heart.av.Aann*system.heart.lv.E[n+2]/(system.solverparams.rho*
        #     system.heart.av.leff*Vs*zs)*x[3]^2+ts*Vs/(2*system.heart.av.Aann*
        #     system.heart.av.leff*system.solverparams.h*zs^3)*((-system.heart.lv.V[n+1]/Vs+
        #     x[2])-abs(system.heart.lv.V[n+1]/Vs-x[2]))+ts*Vs*system.heart.av.Aajnn*
        #     system.heart.av.Ks*system.heart.lv.E[n+2]/(system.solverparams.rho*
        #     system.heart.av.leff*zs)*x[3]^2*(system.heart.lv.V[n+1]/Vs+
        #     system.heart.lv.V0/Vs-2*x[2]));
    end
    if state == "opening"
        J33 = (3*x[3]^2*(system.heart.lv.V[n+1]/Vs-x[2])*ts/system.solverparams.h-
            system.heart.av.Aann*ts/(system.solverparams.rho*system.heart.av.leff*
            system.heart.av.Kvo*Vs*zs^2)*((system.heart.av.zeta[n+1]*x[3]^2-zs*x[3]^3)/
            ((zs^-1-x[3])^2)+(3*zs*x[3]^2-2*system.heart.av.zeta[n+1]*x[3])/(zs^-1-x[3]))+
            2*system.heart.av.Ks*system.heart.lv.E[n+2]*system.heart.av.Aann*ts*Vs/
            (system.solverparams.rho*system.heart.av.leff*zs)*
            (x[2]-system.heart.lv.V0/Vs)*(system.heart.lv.V[n+1]/Vs-x[2])*x[3]-
            3*x[3]^2*ts/Vs*system.branches.Q[1][n+1,1]);
    else
        J33 = (3*x[3]^2*(system.heart.lv.V[n+1]/Vs-x[2])*ts/system.solverparams.h-
            system.heart.av.Aann*ts/(system.solverparams.rho*system.heart.av.leff*
            system.heart.av.Kvc*Vs*zs^2)*(2*x[3]*zs-system.heart.av.zeta[n+1])+
            2*system.heart.av.Ks*system.heart.lv.E[n+2]*system.heart.av.Aann*ts*Vs/
            (system.solverparams.rho*system.heart.av.leff*zs)*
            (x[2]-system.heart.lv.V0/Vs)*(system.heart.lv.V[n+1]/Vs-x[2])*x[3]-
            3*x[3]^2*ts/Vs*system.branches.Q[1][n+1,1]);
    end
    # J33 = 3*x[3]^2*(system.heart.lv.V[n+1]/Vs-x[2])/system.solverparams.h*ts-
    #     ts*system.solverparams.h*system.heart.av.Aann/(system.solverparams.rho*
    #     system.heart.av.leff*Vs)*(2*x[3]/zs*(system.heart.lv.E[n+2]*Vs*
    #     (x[2]-system.heart.lv.V0/Vs)-system.branches.beta[1][end]*
    #     ((2*system.solverparams.rho/system.branches.beta[1][end])*vs^2*
    #     ((x[1]-system.branches.W2root/vs)/8+system.branches.c0[1][end])^2-
    #     system.branches.A0[1][end]))-
    #     system.heart.av.Ks*system.heart.lv.E[n+2]*Vs^2/
    #     (zs*system.solverparams.h)*2*x[3]*(x[2]-system.heart.lv.V0/Vs)*
    #     (system.heart.lv.V[n+1]/Vs-x[2]))-3*x[3]^2*ts/Vs*system.branches.Q[1][n+1,1];
    J3 = [J31 J32 J33];

    J = [J1;J2;J3];
end
