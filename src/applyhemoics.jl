function applyhemoics!(system::CVSystem,old=Dict("a"=>0))
    for i = 1:length(system.branches.ID)
        temp = false;
        for j = 1:length(system.hemo.hID)
            if system.branches.ID[i] == system.hemo.hID[j]
                temp = true;
            end
        end
        push!(system.hemo.injured,temp);
    end
    if haskey(old,"hemo")
        hemo = old["hemo"];
        if hemo["Vloss"] == 0
            for j = 1:length(system.hemo.tID)
                push!(system.hemo.bmin,system.branches.beta[system.hemo.tID[j]][end]);
                push!(system.hemo.bmax,2*system.branches.beta[system.hemo.tID[j]][end]);
                push!(system.hemo.Amax,system.branches.A0[system.hemo.tID[j]][end]);
                push!(system.hemo.Amin,0.01*system.branches.A0[system.hemo.tID[j]][end]);
            end
        elseif hemo["Vloss"] > 0
            for j = 1:length(system.hemo.tID)
                push!(system.hemo.bmin,hemo["bmin"][j]);
                push!(system.hemo.bmax,hemo["bmax"][j]);
                push!(system.hemo.Amin,hemo["Amin"][j]);
                push!(system.hemo.Amax,hemo["Amax"][j]);
            end
            system.hemo.Vloss = hemo["Vloss"];
            system.hemo.Vlossinit = system.hemo.Vloss;
        end
    else
        for j = 1:length(system.hemo.tID)
            push!(system.hemo.bmin,system.branches.beta[system.hemo.tID[j]][end]);
            push!(system.hemo.bmax,2*system.branches.beta[system.hemo.tID[j]][end]);
            push!(system.hemo.Amax,system.branches.A0[system.hemo.tID[j]][end]);
            push!(system.hemo.Amin,0.01*system.branches.A0[system.hemo.tID[j]][end]);
        end
    end
end
