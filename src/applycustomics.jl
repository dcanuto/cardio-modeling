function applycustomics!(system::CVSystem)
    # custom peripheral volumes for better initial blood distribution
    system.branches.term[13].P[1,2] = 65*mmHgToPa;
    system.branches.term[13].V[1,2] = (system.branches.term[13].P[1,2]*
        system.branches.term[13].C[2] + system.branches.term[13].V0[2]);
    system.branches.term[13].P[1,3] = 32*mmHgToPa;
    system.branches.term[13].V[1,3] = (system.branches.term[13].P[1,3]*
        system.branches.term[13].C[3] + system.branches.term[13].V0[3]);
    system.branches.term[14].P[1,2] = 73*mmHgToPa;
    system.branches.term[14].V[1,2] = (system.branches.term[14].P[1,2]*
        system.branches.term[14].C[2] + system.branches.term[14].V0[2]);
    system.branches.term[14].P[1,3] = 36*mmHgToPa;
    system.branches.term[14].V[1,3] = (system.branches.term[14].P[1,3]*
        system.branches.term[14].C[3] + system.branches.term[14].V0[3]);
    system.branches.term[16].P[1,2] = 78*mmHgToPa;
    system.branches.term[16].V[1,2] = (system.branches.term[16].P[1,2]*
        system.branches.term[16].C[2] + system.branches.term[16].V0[2]);
    system.branches.term[16].P[1,3] = 38.5*mmHgToPa;
    system.branches.term[16].V[1,3] = (system.branches.term[16].P[1,3]*
        system.branches.term[16].C[3] + system.branches.term[16].V0[3]);
    system.branches.term[18].P[1,2] = 78*mmHgToPa;
    system.branches.term[18].V[1,2] = (system.branches.term[18].P[1,2]*
        system.branches.term[18].C[2] + system.branches.term[18].V0[2]);
    system.branches.term[18].P[1,3] = 38.5*mmHgToPa;
    system.branches.term[18].V[1,3] = (system.branches.term[18].P[1,3]*
        system.branches.term[18].C[3] + system.branches.term[18].V0[3]);
    system.branches.term[21].P[1,2] = 60*mmHgToPa;
    system.branches.term[21].V[1,2] = (system.branches.term[21].P[1,2]*
        system.branches.term[21].C[2] + system.branches.term[21].V0[2]);
    system.branches.term[21].P[1,3] = 30*mmHgToPa;
    system.branches.term[21].V[1,3] = (system.branches.term[21].P[1,3]*
        system.branches.term[21].C[3] + system.branches.term[21].V0[3]);
    system.branches.term[23].P[1,2] = 60*mmHgToPa;
    system.branches.term[23].V[1,2] = (system.branches.term[23].P[1,2]*
        system.branches.term[23].C[2] + system.branches.term[23].V0[2]);
    system.branches.term[23].P[1,3] = 30*mmHgToPa;
    system.branches.term[23].V[1,3] = (system.branches.term[23].P[1,3]*
        system.branches.term[23].C[3] + system.branches.term[23].V0[3]);
    system.branches.term[27].P[1,2] = 61*mmHgToPa;
    system.branches.term[27].V[1,2] = (system.branches.term[27].P[1,2]*
        system.branches.term[27].C[2] + system.branches.term[27].V0[2]);
    system.branches.term[27].P[1,3] = 30.5*mmHgToPa;
    system.branches.term[27].V[1,3] = (system.branches.term[27].P[1,3]*
        system.branches.term[27].C[3] + system.branches.term[27].V0[3]);
    system.branches.term[33].P[1,2] = 42*mmHgToPa;
    system.branches.term[33].V[1,2] = (system.branches.term[33].P[1,2]*
        system.branches.term[33].C[2] + system.branches.term[33].V0[2]);
    system.branches.term[33].P[1,3] = 21*mmHgToPa;
    system.branches.term[33].V[1,3] = (system.branches.term[33].P[1,3]*
        system.branches.term[33].C[3] + system.branches.term[33].V0[3]);
    system.branches.term[35].P[1,2] = 16.4*mmHgToPa;
    system.branches.term[35].V[1,2] = (system.branches.term[35].P[1,2]*
        system.branches.term[35].C[2] + system.branches.term[35].V0[2]);
    system.branches.term[35].P[1,3] = 9.7*mmHgToPa;
    system.branches.term[35].V[1,3] = (system.branches.term[35].P[1,3]*
        system.branches.term[35].C[3] + system.branches.term[35].V0[3]);
    system.branches.term[37].P[1,2] = 60*mmHgToPa;
    system.branches.term[37].V[1,2] = (system.branches.term[37].P[1,2]*
        system.branches.term[37].C[2] + system.branches.term[37].V0[2]);
    system.branches.term[37].P[1,3] = 29*mmHgToPa;
    system.branches.term[37].V[1,3] = (system.branches.term[37].P[1,3]*
        system.branches.term[37].C[3] + system.branches.term[37].V0[3]);
    system.branches.term[39].P[1,2] = 62*mmHgToPa;
    system.branches.term[39].V[1,2] = (system.branches.term[39].P[1,2]*
        system.branches.term[39].C[2] + system.branches.term[39].V0[2]);
    system.branches.term[39].P[1,3] = 30.5*mmHgToPa;
    system.branches.term[39].V[1,3] = (system.branches.term[39].P[1,3]*
        system.branches.term[39].C[3] + system.branches.term[39].V0[3]);
    system.branches.term[47].P[1,2] = 16.3*mmHgToPa;
    system.branches.term[47].V[1,2] = (system.branches.term[47].P[1,2]*
        system.branches.term[47].C[2] + system.branches.term[47].V0[2]);
    system.branches.term[47].P[1,3] = 9.7*mmHgToPa;
    system.branches.term[47].V[1,3] = (system.branches.term[47].P[1,3]*
        system.branches.term[47].C[3] + system.branches.term[47].V0[3]);
    system.branches.term[48].P[1,2] = 42*mmHgToPa;
    system.branches.term[48].V[1,2] = (system.branches.term[48].P[1,2]*
        system.branches.term[48].C[2] + system.branches.term[48].V0[2]);
    system.branches.term[48].P[1,3] = 20.8*mmHgToPa;
    system.branches.term[48].V[1,3] = (system.branches.term[48].P[1,3]*
        system.branches.term[48].C[3] + system.branches.term[48].V0[3]);
    system.branches.term[90].P[1,2] = 16.6*mmHgToPa;
    system.branches.term[90].V[1,2] = (system.branches.term[90].P[1,2]*
        system.branches.term[90].C[2] + system.branches.term[90].V0[2]);
    system.branches.term[90].P[1,3] = 9.9*mmHgToPa;
    system.branches.term[90].V[1,3] = (system.branches.term[90].P[1,3]*
        system.branches.term[90].C[3] + system.branches.term[90].V0[3]);
    system.branches.term[54].P[1,2] = 16.4*mmHgToPa;
    system.branches.term[54].V[1,2] = (system.branches.term[54].P[1,2]*
        system.branches.term[54].C[2] + system.branches.term[54].V0[2]);
    system.branches.term[54].P[1,3] = 9.75*mmHgToPa;
    system.branches.term[54].V[1,3] = (system.branches.term[54].P[1,3]*
        system.branches.term[54].C[3] + system.branches.term[54].V0[3]);
    system.branches.term[55].P[1,2] = 16.4*mmHgToPa;
    system.branches.term[55].V[1,2] = (system.branches.term[55].P[1,2]*
        system.branches.term[55].C[2] + system.branches.term[55].V0[2]);
    system.branches.term[55].P[1,3] = 9.75*mmHgToPa;
    system.branches.term[55].V[1,3] = (system.branches.term[55].P[1,3]*
        system.branches.term[55].C[3] + system.branches.term[55].V0[3]);
    system.branches.term[57].P[1,2] = 60*mmHgToPa;
    system.branches.term[57].V[1,2] = (system.branches.term[57].P[1,2]*
        system.branches.term[57].C[2] + system.branches.term[57].V0[2]);
    system.branches.term[57].P[1,3] = 29.2*mmHgToPa;
    system.branches.term[57].V[1,3] = (system.branches.term[57].P[1,3]*
        system.branches.term[57].C[3] + system.branches.term[57].V0[3]);
    system.branches.term[58].P[1,2] = 57*mmHgToPa;
    system.branches.term[58].V[1,2] = (system.branches.term[58].P[1,2]*
        system.branches.term[58].C[2] + system.branches.term[58].V0[2]);
    system.branches.term[58].P[1,3] = 27.7*mmHgToPa;
    system.branches.term[58].V[1,3] = (system.branches.term[58].P[1,3]*
        system.branches.term[58].C[3] + system.branches.term[58].V0[3]);
    system.branches.term[67].P[1,2] = 44*mmHgToPa;
    system.branches.term[67].V[1,2] = (system.branches.term[67].P[1,2]*
        system.branches.term[67].C[2] + system.branches.term[67].V0[2]);
    system.branches.term[67].P[1,3] = 21.7*mmHgToPa;
    system.branches.term[67].V[1,3] = (system.branches.term[67].P[1,3]*
        system.branches.term[67].C[3] + system.branches.term[67].V0[3]);
    system.branches.term[68].P[1,2] = 18.2*mmHgToPa;
    system.branches.term[68].V[1,2] = (system.branches.term[68].P[1,2]*
        system.branches.term[68].C[2] + system.branches.term[68].V0[2]);
    system.branches.term[68].P[1,3] = 10.5*mmHgToPa;
    system.branches.term[68].V[1,3] = (system.branches.term[68].P[1,3]*
        system.branches.term[68].C[3] + system.branches.term[68].V0[3]);
    system.branches.term[70].P[1,2] = 54*mmHgToPa;
    system.branches.term[70].V[1,2] = (system.branches.term[70].P[1,2]*
        system.branches.term[70].C[2] + system.branches.term[70].V0[2]);
    system.branches.term[70].P[1,3] = 26.5*mmHgToPa;
    system.branches.term[70].V[1,3] = (system.branches.term[70].P[1,3]*
        system.branches.term[70].C[3] + system.branches.term[70].V0[3]);
    system.branches.term[73].P[1,2] = 58.8*mmHgToPa;
    system.branches.term[73].V[1,2] = (system.branches.term[73].P[1,2]*
        system.branches.term[73].C[2] + system.branches.term[73].V0[2]);
    system.branches.term[73].P[1,3] = 28.8*mmHgToPa;
    system.branches.term[73].V[1,3] = (system.branches.term[73].P[1,3]*
        system.branches.term[73].C[3] + system.branches.term[73].V0[3]);
    system.branches.term[83].P[1,2] = 17.9*mmHgToPa;
    system.branches.term[83].V[1,2] = (system.branches.term[83].P[1,2]*
        system.branches.term[83].C[2] + system.branches.term[83].V0[2]);
    system.branches.term[83].P[1,3] = 10.3*mmHgToPa;
    system.branches.term[83].V[1,3] = (system.branches.term[83].P[1,3]*
        system.branches.term[83].C[3] + system.branches.term[83].V0[3]);
    system.branches.term[84].P[1,2] = 54*mmHgToPa;
    system.branches.term[84].V[1,2] = (system.branches.term[84].P[1,2]*
        system.branches.term[84].C[2] + system.branches.term[84].V0[2]);
    system.branches.term[84].P[1,3] = 26.4*mmHgToPa;
    system.branches.term[84].V[1,3] = (system.branches.term[84].P[1,3]*
        system.branches.term[84].C[3] + system.branches.term[84].V0[3]);
    system.branches.term[85].P[1,2] = 43.3*mmHgToPa;
    system.branches.term[85].V[1,2] = (system.branches.term[85].P[1,2]*
        system.branches.term[85].C[2] + system.branches.term[85].V0[2]);
    system.branches.term[85].P[1,3] = 21.5*mmHgToPa;
    system.branches.term[85].V[1,3] = (system.branches.term[85].P[1,3]*
        system.branches.term[85].C[3] + system.branches.term[85].V0[3]);
    system.branches.term[91].P[1,2] = 16*mmHgToPa;
    system.branches.term[91].V[1,2] = (system.branches.term[91].P[1,2]*
        system.branches.term[91].C[2] + system.branches.term[91].V0[2]);
    system.branches.term[91].P[1,3] = 9.52*mmHgToPa;
    system.branches.term[91].V[1,3] = (system.branches.term[91].P[1,3]*
        system.branches.term[91].C[3] + system.branches.term[91].V0[3]);
end
