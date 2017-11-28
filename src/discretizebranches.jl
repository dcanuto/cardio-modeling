function discretizebranches!(system::CVSystem)
    h = [];

    # branch grid spacing
    for i in 1:length(system.branches.ID)
        push!(system.branches.k,system.branches.lengthincm[i]*cmTom/
            (system.solverparams.JL-1));
        push!(h,system.solverparams.CFL*system.branches.k[i]/
            system.branches.c0[i][end]);
    end

    # time step guaranteed to satisfy CFL for all branches
    system.solverparams.h = minimum(h);
    system.solverparams.numsteps = ceil(system.heart.activation.th[1]/
        system.solverparams.h);

    # allocate space for 1D domain solution variables
    for i in 1:length(system.branches.ID)
        push!(system.branches.A,
            zeros(system.solverparams.numsteps+1,system.solverparams.JL));
        push!(system.branches.Q,
            zeros(system.solverparams.numsteps+1,system.solverparams.JL));
        push!(system.branches.P,
            zeros(system.solverparams.numsteps+1,system.solverparams.JL));
        push!(system.branches.Fp,
            zeros(2*system.solverparams.JL));
        push!(system.branches.Fbarforward,
            zeros(2*system.solverparams.JL-4));
        push!(system.branches.Fbarbackward,
            zeros(2*system.solverparams.JL-4));
        push!(system.branches.Abackward,
            zeros(system.solverparams.JL-2));
        push!(system.branches.Aforward,
            zeros(system.solverparams.JL-2));
        push!(system.branches.Qbackward,
            zeros(system.solverparams.JL-2));
        push!(system.branches.Qforward,
            zeros(system.solverparams.JL-2));
        push!(system.branches.W1end,0.);
        push!(system.branches.W1,0.);
        push!(system.branches.W2,0.);
    end

    # discretize time for first cardiac cycle
    system.t = system.solverparams.h*[0:1:size(system.branches.A[1],1)-1;];

    return system
end
