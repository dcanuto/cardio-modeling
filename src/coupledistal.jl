function coupledistal!(system::CVSystem,n::Int64)
    # W1 at next time step
    endinvariants!(system,n);

    for i = 1:length(system.branches.ID)
        if isempty(system.branches.children[i])
            newtondist!(system,n,i);
        end
    end
end
