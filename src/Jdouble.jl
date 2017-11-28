function Jdouble(x::Vector{Float64},system::CVSystem,parentID::Int64,children::Vector{Int64})
    J11 = 1;
    J21 = system.solverparams.rho*x[1]/(x[2]^2);
    J31 = J21;
    J41 = 1/x[2];
    J51 = 0;
    J61 = 0;
    J1 = [J11;J21;J31;J41;J51;J61];

    J12 = 0;
    J22 = (0.5*system.branches.beta[parentID][end]*x[2]^-0.5-
        system.solverparams.rho*x[1]^2*x[2]^-3);
    J32 = J22;
    J42 = -x[1]*x[2]^-2 + sqrt(0.5*system.branches.beta[parentID][end]/
        system.solverparams.rho)*x[2]^-0.75;
    J52 = 0;
    J62 = 0;
    J2 = [J12;J22;J32;J42;J52;J62];

    J13 = -1;
    J23 = -system.solverparams.rho*x[3]/(x[4]^2);
    J33 = 0;
    J43 = 0;
    J53 = 1/x[4];
    J63 = 0;
    J3 = [J13;J23;J33;J43;J53;J63];

    J14 = 0;
    J24 = (-0.5*system.branches.beta[children[1]][end]*x[4]^-0.5+
        system.solverparams.rho*x[3]^2*x[4]^-3);
    J34 = 0;
    J44 = 0;
    J54 = -x[3]*x[4]^-2 - sqrt(0.5*system.branches.beta[children[1]][end]/
        system.solverparams.rho)*x[4]^-0.75;
    J64 = 0;
    J4 = [J14;J24;J34;J44;J54;J64];

    J15 = -1;
    J25 = 0;
    J35 = -system.solverparams.rho*x[5]/(x[6]^2);
    J45 = 0;
    J55 = 0;
    J65 = 1/x[6];
    J5 = [J15;J25;J35;J45;J55;J65];

    J16 = 0;
    J26 = 0;
    J36 = (-0.5*system.branches.beta[children[2]][end]*x[6]^-0.5+
        system.solverparams.rho*x[5]^2*x[6]^-3);
    J46 = 0;
    J56 = 0;
    J66 = -x[5]*x[6]^-2 - sqrt(0.5*system.branches.beta[children[2]][end]/
        system.solverparams.rho)*x[6]^-0.75;
    J6 = [J16;J26;J36;J46;J56;J66];

    J = [J1 J2 J3 J4 J5 J6];
end
