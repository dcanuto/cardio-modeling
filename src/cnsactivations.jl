function cnsactivations!(system::CVSystem,n::Int64)
    # activation based on Gompertz's sigmoidal function
    push!(system.cns.np,system.cns.a*exp(
        -system.cns.b*exp(-system.cns.c*
        (system.cns.Paverage[system.solverparams.numbeats]/
        system.cns.Ptarget - 1))));
    push!(system.cns.ns,system.cns.a*exp(
        -system.cns.b*exp(system.cns.c*
        (system.cns.Paverage[system.solverparams.numbeats]/
        system.cns.Ptarget - 1))));
end
