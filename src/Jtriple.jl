function Jtriple(x::Vector{Float64},system::CVSystem,parentID::Int64,children::Vector{Int64})
    J11 = 1;
    J21 = system.solverparams.rho*x[1]/(x[2]^2);
    J31 = J21;
    J41 = J21;
    J51 = 1/x[2];
    J61 = 0;
    J71 = 0;
    J81 = 0;
    J1 = [J11;J21;J31;J41;J51;J61;J71;J81];

    J12 = 0;
    J22 = (0.5*system.branches.beta[parentID][end]*x[2]^-0.5-
        system.solverparams.rho*x[1]^2*x[2]^-3);
    J32 = J22;
    J42 = J22;
    J52 = -x[1]*x[2]^-2 + sqrt(0.5*system.branches.beta[parentID][end]/
        system.solverparams.rho)*x[2]^-0.75;
    J62 = 0;
    J72 = 0;
    J82 = 0;
    J2 = [J12;J22;J32;J42;J52;J62;J72;J82];

    J13 = -1;
    J23 = -system.solverparams.rho*x[3]/(x[4]^2);
    J33 = 0;
    J43 = 0;
    J53 = 0;
    J63 = 1/x[4];
    J73 = 0;
    J83 = 0;
    J3 = [J13;J23;J33;J43;J53;J63;J73;J83];

    J14 = 0;
    J24 = (-0.5*system.branches.beta[children[1]][end]*x[4]^-0.5+
        system.solverparams.rho*x[3]^2*x[4]^-3);
    J34 = 0;
    J44 = 0;
    J54 = 0;
    J64 = -x[3]*x[4]^-2 - sqrt(0.5*system.branches.beta[children[1]][end]/
        system.solverparams.rho)*x[4]^-0.75;
    J74 = 0;
    J84 = 0;
    J4 = [J14;J24;J34;J44;J54;J64;J74;J84];

    J15 = -1;
    J25 = 0;
    J35 = -system.solverparams.rho*x[5]/(x[6]^2);
    J45 = 0;
    J55 = 0;
    J65 = 0;
    J75 = 1/x[6];
    J85 = 0;
    J5 = [J15;J25;J35;J45;J55;J65;J75;J85];

    J16 = 0;
    J26 = 0;
    J36 = (-0.5*system.branches.beta[children[2]][end]*x[6]^-0.5+
        system.solverparams.rho*x[5]^2*x[6]^-3);
    J46 = 0;
    J56 = 0;
    J66 = 0;
    J76 = -x[5]*x[6]^-2 - sqrt(0.5*system.branches.beta[children[2]][end]/
        system.solverparams.rho)*x[6]^-0.75;
    J86 = 0;
    J6 = [J16;J26;J36;J46;J56;J66;J76;J86];

    J17 = -1;
    J27 = 0;
    J37 = 0;
    J47 = -system.solverparams.rho*x[7]/(x[8]^2);
    J57 = 0;
    J67 = 0;
    J77 = 0;
    J87 = 1/x[8];
    J7 = [J17;J27;J37;J47;J57;J67;J77;J87];

    J18 = 0;
    J28 = 0;
    J38 = 0;
    J48 = (-0.5*system.branches.beta[children[3]][end]*x[8]^-0.5+
        system.solverparams.rho*x[7]^2*x[8]^-3);
    J58 = 0;
    J68 = 0;
    J78 = 0;
    J88 = -x[7]*x[8]^-2 - sqrt(0.5*system.branches.beta[children[3]][end]/
        system.solverparams.rho)*x[8]^-0.75;
    J8 = [J18;J28;J38;J48;J58;J68;J78;J88];

    J = [J1 J2 J3 J4 J5 J6 J7 J8];
end
