function fav(x::Vector{Float64},system::CVSystem,n::Int64,state::String)
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

    if length(x) == 3
        f1 = ((system.heart.lv.V[n+1]/Vs - x[2])*ts/system.solverparams.h-ts*vs^5/Vs*
            (2*system.solverparams.rho/system.branches.beta[1][end])^2*
            ((x[1] - system.branches.W2root/vs)/8 + system.branches.c0[1][end]/vs)^4*
            ((x[1] + system.branches.W2root/vs)/2));
        if state == "opening"
            f2 = (x[3] - system.solverparams.h*((1/zs-x[3])*system.heart.av.Kvo*
                (system.heart.lv.E[n+2]*Vs*(x[2]-system.heart.lv.V0/Vs)-
                system.branches.beta[1][end]*((2*system.solverparams.rho/
                system.branches.beta[1][end])*vs^2*((x[1]-system.branches.W2root/vs)/8+
                system.branches.c0[1][end]/vs)^2-sqrt(system.branches.A0[1][end]))))-
                system.heart.av.zeta[n+1]/zs);
            f3 = (x[3]*(system.heart.lv.V[n+1]/Vs-x[2])*ts/system.solverparams.h-
                system.heart.av.Aann*ts/(system.solverparams.rho*system.heart.av.leff*Vs)*
                ((zs*x[3]^3-system.heart.av.zeta[n+1]*(x[3])^2)/(system.heart.av.Kvo*
                (1/zs-x[3]))-1/zs*system.solverparams.rho*Vs^2/(2*system.heart.av.Aann^2*
                system.solverparams.h)*(system.heart.lv.V[n+1]/Vs-x[2])*
                abs(system.heart.lv.V[n+1]/Vs-x[2])-zs*x[3]^2*system.heart.av.Ks*
                system.heart.lv.E[n+2]*Vs^2*(x[2]-system.heart.lv.V0/Vs)*
                (system.heart.lv.V[n+1]/Vs-x[2]))-
                ts/Vs*x[3]*system.branches.Q[1][n+1,1]);
        else
            f2 = (x[3] - system.solverparams.h*((x[3])*system.heart.av.Kvc*
                (system.heart.lv.E[n+2]*Vs*(x[2]-system.heart.lv.V0/Vs)-
                system.branches.beta[1][end]*((2*system.solverparams.rho/
                system.branches.beta[1][end])*vs^2*((x[1]-system.branches.W2root/vs)/8+
                system.branches.c0[1][end]/vs)^2-sqrt(system.branches.A0[1][end]))))-
                system.heart.av.zeta[n+1]/zs);
            f3 = (x[3]*(system.heart.lv.V[n+1]/Vs-x[2])*ts/system.solverparams.h-
                system.heart.av.Aann*ts/(system.solverparams.rho*system.heart.av.leff*Vs)*
                ((zs*x[3]^2-system.heart.av.zeta[n+1]*x[3])/system.heart.av.Kvc
                -1/zs*system.solverparams.rho*Vs^2/(2*system.heart.av.Aann^2*
                system.solverparams.h)*(system.heart.lv.V[n+1]/Vs-x[2])*
                abs(system.heart.lv.V[n+1]/Vs-x[2])-zs*x[3]^2*system.heart.av.Ks*
                system.heart.lv.E[n+2]*Vs^2*(x[2]-system.heart.lv.V0/Vs)*
                (system.heart.lv.V[n+1]/Vs-x[2]))-
                ts/Vs*x[3]*system.branches.Q[1][n+1,1]);
        end
        f = [f1,f2,f3];
    else
        f1 = ((system.heart.lv.V[n+1]/Vs - x[2])*ts/system.solverparams.h-ts*vs^5/Vs*
            (2*system.solverparams.rho/system.branches.beta[1][end])^2*
            ((x[1] - system.branches.W2root/vs)/8 + system.branches.c0[1][end]/vs)^4*
            ((x[1] + system.branches.W2root/vs)/2));
        f3 = ((system.heart.lv.V[n+1]/Vs-x[2])*ts/system.solverparams.h-
            system.heart.av.Aann*ts/(system.solverparams.rho*system.heart.av.leff*Vs)*
            (-1/zs*system.solverparams.rho*Vs^2/(2*system.heart.av.Aann^2*
            system.solverparams.h)*(system.heart.lv.V[n+1]/Vs-x[2])*
            abs(system.heart.lv.V[n+1]/Vs-x[2])-zs*system.heart.av.Ks*
            system.heart.lv.E[n+2]*Vs^2*(x[2]-system.heart.lv.V0/Vs)*
            (system.heart.lv.V[n+1]/Vs-x[2]))-
            ts/Vs*system.branches.Q[1][n+1,1]);
        f = [f1,f3];
    end
end
