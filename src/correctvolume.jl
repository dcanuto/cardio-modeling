function correctvolume!(system::CVSystem,n::Int64)
    # get current volume
    updatevolumes!(system,n);

    # find error in volume
    error = system.finalvolume - system.initialvolume;

    # divide error proportionally to accumulation regions
    lungerror1 = 0.05*error;
    lungerror2 = 0.25*error;
    termerror = 0.5*error;
    lverror = 0.1*error;
    rverror = 0.1*error;

    # subtract erroneous volume to maintain mass conservation
    system.lungs.Vv[n+2,1] -= lungerror1;
    system.lungs.Vv[n+2,2] -= lungerror2;
    system.heart.lv.V[n+2] -= lverror;
    system.heart.rv.V[n+2] -= rverror;

    numterms = 0;
    for i = 1:length(system.branches.ID)
        if isempty(system.branches.children[i])
            numterms+=1;
        end
    end
    termerror = termerror/numterms;
    for i = 1:length(system.branches.ID)
        if isempty(system.branches.children[i])
            system.branches.term[i].V[n+2,5] -= termerror;
        end
    end
end
