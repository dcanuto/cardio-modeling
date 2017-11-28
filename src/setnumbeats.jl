function setnumbeats!(system::CVSystem,n::Int64)
    # total time needed for all past/current cardiac cycles
    ttotal = sum(system.heart.activation.th);

    # increment number of cycles if new cycle is starting
    if system.t[n+2] >= ttotal
        system.solverparams.numbeats+=1;
    end
end
