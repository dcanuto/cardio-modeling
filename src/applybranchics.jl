function applybranchics!(system::CVSystem,old=Dict("a"=>0),restart="no")
    if restart == "no"
        Pstart = 70*mmHgToPa;
        for i in 1:length(system.branches.ID)
            system.branches.A[i][1,:] = (Pstart/system.branches.beta[i][end] +
                sqrt(system.branches.A0[i][end]))^2;
            system.branches.Q[i][1,:] = 0;
            system.branches.P[i][1,:] = Pstart/mmHgToPa;
        end
    elseif restart == "yes"
        for i in 1:length(system.branches.ID)
            system.branches.A[i][1,:] = old["A"][i][end,:];
            system.branches.Q[i][1,:] = old["Q"][i][end,:];
            system.branches.P[i][1,:] = old["P"][i][end,:];
        end
    end
end
