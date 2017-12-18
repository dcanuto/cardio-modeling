function applylungics!(system::CVSystem,old=Dict("a"=>0),restart="no")
    if restart == "no"
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
    elseif restart == "yes"
        system.lungs.Pp[1] = old["Pp"][end];
        system.lungs.Pa[1,1] = old["Pa"][end,1];
        system.lungs.Va[1,1] = old["Va"][end,1];
        system.lungs.Qa[1,1] = old["Qa"][end,1];
        system.lungs.Pa[1,2] = old["Pa"][end,2];
        system.lungs.Va[1,2] = old["Va"][end,2];
        system.lungs.Qa[1,2] = old["Qa"][end,2];
        system.lungs.Pa[1,3] = old["Pa"][end,3];
        system.lungs.Va[1,3] = old["Va"][end,3];
        system.lungs.Qa[1,3] = old["Qa"][end,3];
        system.lungs.Pv[1,1] = old["Pv"][end,1];
        system.lungs.Vv[1,1] = old["Vv"][end,1];
        system.lungs.Qv[1,1] = old["Qv"][end,1];
        system.lungs.Pv[1,2] = old["Pv"][end,2];
        system.lungs.Vv[1,2] = old["Vv"][end,2];
        system.lungs.Qv[1,2] = old["Qv"][end,2];
    end
end
