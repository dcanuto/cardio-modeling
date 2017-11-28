function discretizelungs!(system::CVSystem)
    append!(system.lungs.Pp,zeros(system.solverparams.numsteps+1))

    system.lungs.Pa = zeros(system.solverparams.numsteps+1,3);
    system.lungs.Va = zeros(system.solverparams.numsteps+1,3);
    system.lungs.Qa = zeros(system.solverparams.numsteps+1,3);
    system.lungs.Pv = zeros(system.solverparams.numsteps+1,2);
    system.lungs.Vv = zeros(system.solverparams.numsteps+1,2);
    system.lungs.Qv = zeros(system.solverparams.numsteps+1,2);
end
