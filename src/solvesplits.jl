function solvesplits!(system::CVSystem,n::Int64)
    for i = 1:length(system.branches.ID)
        if ~isempty(system.branches.children[i])
            newton!(system,n,i);
        end
    end
end
