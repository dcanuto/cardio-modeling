function applylungics!(system::CVSystem)
    system.lungs.Pp[1] = 8*mmHgToPa;
    system.lungs.Pa[1,1] = 7.5*mmHgToPa;
    system.lungs.Va[1,1] = (system.lungs.Pa[1,1]*system.lungs.Ca[1] +
        system.lungs.V0[1]);
    system.lungs.Pa[1,2] = 7*mmHgToPa;
    system.lungs.Va[1,2] = (system.lungs.Pa[1,2]*system.lungs.Ca[2] +
        system.lungs.V0[2]);
    system.lungs.Pa[1,3] = 6.5*mmHgToPa;
    system.lungs.Va[1,3] = (system.lungs.Pa[1,3]*system.lungs.Ca[3] +
        system.lungs.V0[3]);
    system.lungs.Pv[1,1] = 4.5*mmHgToPa;
    system.lungs.Vv[1,1] = (system.lungs.Pv[1,1]*system.lungs.Cv[1] +
        system.lungs.V0[4]);
    system.lungs.Pv[1,2] = 4.5*mmHgToPa;
    system.lungs.Vv[1,2] = (system.lungs.Pv[1,2]*system.lungs.Cv[2] +
        system.lungs.V0[5]);    
end
