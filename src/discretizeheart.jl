function discretizeheart!(system::CVSystem)
    append!(system.heart.lv.V,zeros(system.solverparams.numsteps+1))
    append!(system.heart.lv.P,zeros(system.solverparams.numsteps+1))
    append!(system.heart.lv.E,zeros(system.solverparams.numsteps+1))

    append!(system.heart.rv.V,zeros(system.solverparams.numsteps+1))
    append!(system.heart.rv.P,zeros(system.solverparams.numsteps+1))
    append!(system.heart.rv.Q,zeros(system.solverparams.numsteps+1))
    append!(system.heart.rv.E,zeros(system.solverparams.numsteps+1))

    append!(system.heart.la.V,zeros(system.solverparams.numsteps+1))
    append!(system.heart.la.P,zeros(system.solverparams.numsteps+1))
    append!(system.heart.la.Q,zeros(system.solverparams.numsteps+1))

    append!(system.heart.ra.V,zeros(system.solverparams.numsteps+1))
    append!(system.heart.ra.P,zeros(system.solverparams.numsteps+1))
    append!(system.heart.ra.Q,zeros(system.solverparams.numsteps+1))
end
