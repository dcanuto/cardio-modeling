function applycnsics!(system::CVSystem)
    system.cns.H[1] = 1/system.heart.activation.th[1]*minTos;
    system.cns.Emaxlv[1] = system.heart.lv.Emax[1];
    system.cns.Emaxrv[1] = system.heart.rv.Emax[1];

    system.cns.R2L[1] = system.branches.term[13].R[2];
    system.cns.R3L[1] = system.branches.term[13].R[3];
    system.cns.R2U[1] = system.branches.term[50].R[2];
    system.cns.R3U[1] = system.branches.term[50].R[3];
    system.cns.C4L[1] = system.branches.term[13].C[4];
    system.cns.C5L[1] = system.branches.term[13].C[5];
    system.cns.C4U[1] = system.branches.term[50].C[4];
    system.cns.C5U[1] = system.branches.term[50].C[5];
    system.cns.V4L[1] = system.branches.term[13].V0[4];
    system.cns.V5L[1] = system.branches.term[13].V0[5];
    system.cns.V4U[1] = system.branches.term[50].V0[4];
    system.cns.V5U[1] = system.branches.term[50].V0[5];

    system.cns.alphaH = 2.1563;
    system.cns.betaH = 0.5313;
    system.cns.gammaH = 0.8438;
    system.cns.gammalv = 0.8*system.cns.Emaxlv[1];
    system.cns.alphalv = 1.6*system.cns.Emaxlv[1] - system.cns.gammalv;
    system.cns.gammarv = 0.8*system.cns.Emaxrv[1];
    system.cns.alpharv = 1.6*system.cns.Emaxrv[1] - system.cns.gammarv;
    system.cns.gammaR2L = 0.6*system.cns.R2L[1];
    system.cns.alphaR2L = 2.2*system.cns.R2L[1] - system.cns.gammaR2L;
    system.cns.gammaR3L = 0.6*system.cns.R3L[1];
    system.cns.alphaR3L = 2.2*system.cns.R3L[1] - system.cns.gammaR3L;
    system.cns.gammaR2U = 0.6*system.cns.R2U[1];
    system.cns.alphaR2U = 2.2*system.cns.R2U[1] - system.cns.gammaR2U;
    system.cns.gammaR3U = 0.6*system.cns.R3U[1];
    system.cns.alphaR3U = 2.2*system.cns.R3U[1] - system.cns.gammaR3U;
    system.cns.gammaC4L = 1.1*system.cns.C4L[1];
    system.cns.alphaC4L = -0.7*system.cns.C4L[1] + system.cns.gammaC4L;
    system.cns.gammaC5L = 1.1*system.cns.C5L[1];
    system.cns.alphaC5L = -0.7*system.cns.C5L[1] + system.cns.gammaC5L;
    system.cns.gammaC4U = 1.1*system.cns.C4U[1];
    system.cns.alphaC4U = -0.7*system.cns.C4U[1] + system.cns.gammaC4U;
    system.cns.gammaC5U = 1.1*system.cns.C5U[1];
    system.cns.alphaC5U = -0.7*system.cns.C5U[1] + system.cns.gammaC5U;
    system.cns.gammaV4L = 1.05*system.cns.V4L[1];
    system.cns.alphaV4L = -0.85*system.cns.V4L[1] + system.cns.gammaV4L;
    system.cns.gammaV5L = 1.05*system.cns.V5L[1];
    system.cns.alphaV5L = -0.85*system.cns.V5L[1] + system.cns.gammaV5L;
    system.cns.gammaV4U = 1.05*system.cns.V4U[1];
    system.cns.alphaV4U = -0.85*system.cns.V4U[1] + system.cns.gammaV4U;
    system.cns.gammaV5U = 1.05*system.cns.V5U[1];
    system.cns.alphaV5U = -0.85*system.cns.V5U[1] + system.cns.gammaV5U;

    system.cns.ns = [0.25];
    system.cns.np = [0.25];
end
