function discretizebranches!(system::CVSystem,old=Dict("a"=>0),restart="no")
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

    if restart == "no"
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
    elseif restart == "yes"
        # determine time shift
        temp = old["heart"];
        temp = temp["activation"];
        system.solverparams.tshift = old["t"][end] - sum(temp["th"][1:end-1]);

        # update total number of time steps
        ntoadd = ceil((system.heart.activation.th[end]-
            system.solverparams.tshift)/system.solverparams.h);
        system.solverparams.numsteps=ntoadd;

        # update discrete times
        ttoadd = [0:1:ntoadd;]*system.solverparams.h + system.solverparams.tshift;
        system.t = ttoadd;

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
    end

    return system
end
