function solvesplits!(system::CVSystem,n::Int64,hemoflag="no")
    if hemoflag == "no"
        for i = 1:length(system.branches.ID)
            if ~isempty(system.branches.children[i])
                newton!(system,n,i);
            end
        end
    elseif hemoflag == "yes"
        for i = 1:length(system.branches.ID)
            if ~isempty(system.branches.children[i]) && system.hemo.injured[i] == false
                newton!(system,n,i);
            elseif ~isempty(system.branches.children[i]) && system.hemo.injured[i] == true
                modelhemo!(system,n,i);
            end
        end
    end
end
