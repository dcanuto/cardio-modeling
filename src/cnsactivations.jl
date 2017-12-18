function cnsactivations!(system::CVSystem,n::Int64)
    # # activation based on Gompertz's sigmoidal function
    # push!(system.cns.np,system.cns.a*exp(
    #     -system.cns.b*exp(-system.cns.c*
    #     (system.cns.Paverage[system.solverparams.numbeats]/
    #     system.cns.Ptarget - 1))));
    # push!(system.cns.ns,system.cns.a*exp(
    #     -system.cns.b*exp(system.cns.c*
    #     (system.cns.Paverage[system.solverparams.numbeats]/
    #     system.cns.Ptarget - 1))));
    # hard-code activations
    push!(system.cns.ns,0.25);
    push!(system.cns.np,0.25);
end
