function Jsingle(x::Vector{Float64},system::CVSystem,parentID::Int64,children::Vector{Int64})
    J11 = 1.;
    J21 = system.solverparams.rho*x[1]/(x[2]^2);
    J31 = 1/x[2];
    J41 = 0;
    J1 = [J11;J21;J31;J41];

    J12 = 0;
    J22 = (0.5*system.branches.beta[parentID][end]*x[2]^-0.5-
        system.solverparams.rho*x[1]^2*x[2]^-3);
    J32 = -x[1]*x[2]^-2 + sqrt(0.5*system.branches.beta[parentID][end]/
        system.solverparams.rho)*x[2]^-0.75;
    J42 = 0;
    J2 = [J12;J22;J32;J42];

    J13 = -1;
    J23 = -system.solverparams.rho*x[3]/(x[4]^2);
    J33 = 0;
    J43 = 1/x[4];
    J3 = [J13;J23;J33;J43];

    J14 = 0;
    J24 = (-0.5*system.branches.beta[children[1]][end]*x[4]^-0.5+
        system.solverparams.rho*x[3]^2*x[4]^-3);
    J34 = 0;
    J44 = -x[3]*x[4]^-2 - sqrt(0.5*system.branches.beta[children[1]][end]/
        system.solverparams.rho)*x[4]^-0.75;
    J4 = [J14;J24;J34;J44];

    J = [J1 J2 J3 J4];
end
