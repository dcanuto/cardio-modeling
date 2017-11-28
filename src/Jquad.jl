function Jquad(x::Vector{Float64},system::CVSystem,parentID::Int64,children::Vector{Int64})
    J11 = 1;
    J21 = system.solverparams.rho*x[1]/(x[2]^2);
    J31 = J21;
    J41 = J21;
    J51 = J21;
    J61 = 1/x[2];
    J71 = 0;
    J81 = 0;
    J91 = 0;
    J101 = 0;
    J1 = [J11;J21;J31;J41;J51;J61;J71;J81;J91;J101];

    J12 = 0;
    J22 = (0.5*system.branches.beta[parentID][end]*x[2]^-0.5-
        system.solverparams.rho*x[1]^2*x[2]^-3);
    J32 = J22;
    J42 = J22;
    J52 = J22;
    J62 = -x[1]*x[2]^-2 + sqrt(0.5*system.branches.beta[parentID][end]/
        system.solverparams.rho)*x[2]^-0.75;
    J72 = 0;
    J82 = 0;
    J92 = 0;
    J102 = 0;
    J2 = [J12;J22;J32;J42;J52;J62;J72;J82;J92;J102];

    J13 = -1;
    J23 = -system.solverparams.rho*x[3]/(x[4]^2);
    J33 = 0;
    J43 = 0;
    J53 = 0;
    J63 = 0;
    J73 = 1/x[4];
    J83 = 0;
    J93 = 0;
    J103 = 0;
    J3 = [J13;J23;J33;J43;J53;J63;J73;J83;J93;J103];

    J14 = 0;
    J24 = (-0.5*system.branches.beta[children[1]][end]*x[4]^-0.5+
        system.solverparams.rho*x[3]^2*x[4]^-3);
    J34 = 0;
    J44 = 0;
    J54 = 0;
    J64 = 0;
    J74 = -x[3]*x[4]^-2 - sqrt(0.5*system.branches.beta[children[1]][end]/
        system.solverparams.rho)*x[4]^-0.75;
    J84 = 0;
    J94 = 0;
    J104 = 0;
    J4 = [J14;J24;J34;J44;J54;J64;J74;J84;J94;J104];

    J15 = -1;
    J25 = 0;
    J35 = -system.solverparams.rho*x[5]/(x[6]^2);
    J45 = 0;
    J55 = 0;
    J65 = 0;
    J75 = 0;
    J85 = 1/x[6];
    J95 = 0;
    J105 = 0;
    J5 = [J15;J25;J35;J45;J55;J65;J75;J85;J95;J105];

    J16 = 0;
    J26 = 0;
    J36 = (-0.5*system.branches.beta[children[2]][end]*x[6]^-0.5+
        system.solverparams.rho*x[5]^2*x[6]^-3);
    J46 = 0;
    J56 = 0;
    J66 = 0;
    J76 = 0;
    J86 = -x[5]*x[6]^-2 - sqrt(0.5*system.branches.beta[children[2]][end]/
        system.solverparams.rho)*x[6]^-0.75;
    J96 = 0;
    J106 = 0;
    J6 = [J16;J26;J36;J46;J56;J66;J76;J86;J96;J106];

    J17 = -1;
    J27 = 0;
    J37 = 0;
    J47 = -system.solverparams.rho*x[7]/(x[8]^2);
    J57 = 0;
    J67 = 0;
    J77 = 0;
    J87 = 0;
    J97 = 1/x[8];
    J107 = 0;
    J7 = [J17;J27;J37;J47;J57;J67;J77;J87;J97;J107];

    J18 = 0;
    J28 = 0;
    J38 = 0;
    J48 = (-0.5*system.branches.beta[children[3]][end]*x[8]^-0.5+
        system.solverparams.rho*x[7]^2*x[8]^-3);
    J58 = 0;
    J68 = 0;
    J78 = 0;
    J88 = 0;
    J98 = -x[7]*x[8]^-2 - sqrt(0.5*system.branches.beta[children[3]][end]/
        system.solverparams.rho)*x[8]^-0.75;
    J108 = 0;
    J8 = [J18;J28;J38;J48;J58;J68;J78;J88;J98;J108];

    J19 = -1;
    J29 = 0;
    J39 = 0;
    J49 = 0;
    J59 = -system.solverparams.rho*x[9]/(x[10]^2);
    J69 = 0;
    J79 = 0;
    J89 = 0;
    J99 = 0;
    J109 = 1/x[10];
    J9 = [J19;J29;J39;J49;J59;J69;J79;J89;J99;J109];

    J110 = 0;
    J210 = 0;
    J310 = 0;
    J410 = 0;
    J510 = (-0.5*system.branches.beta[children[4]][end]*x[10]^-0.5+
        system.solverparams.rho*x[9]^2*x[10]^-3);
    J610 = 0;
    J710 = 0;
    J810 = 0;
    J910 = 0;
    J1010 = -x[9]*x[10]^-2 - sqrt(0.5*system.branches.beta[children[4]][end]/
        system.solverparams.rho)*x[10]^-0.75;
    J10 = [J110;J210;J310;J410;J510;J610;J710;J810;J910;J1010];

    J = [J1 J2 J3 J4 J5 J6 J7 J8 J9 J10];
end
