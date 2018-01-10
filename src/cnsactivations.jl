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
    nst = 0.25;
    r = 1/system.cns.c*log(log(nst^(-1/system.cns.b)))+1;
    ns = exp(-system.cns.b*exp(system.cns.c*(r-1)));
    np = exp(-system.cns.b*exp(-system.cns.c*(r-1)));
    push!(system.cns.ns,ns);
    push!(system.cns.np,np);
end
