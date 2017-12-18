function applyperipheryics!(system::CVSystem)
    numupper = 0;
    for i = 1:length(system.branches.ID)
        if isempty(system.branches.children[i])
            system.branches.term[i].P[1,1] = (
                system.branches.P[i][1,system.solverparams.JL]*mmHgToPa);
            system.branches.term[i].V[1,1] = (
            system.branches.term[i].P[1,1]*system.branches.term[i].C[1] +
            system.branches.term[i].V0[1]);
            system.branches.term[i].P[1,2] = 80*mmHgToPa;
            system.branches.term[i].V[1,2] = (
            system.branches.term[i].P[1,2]*system.branches.term[i].C[2] +
            system.branches.term[i].V0[2]);
            system.branches.term[i].P[1,3] = 40*mmHgToPa;
            system.branches.term[i].V[1,3] = (
            system.branches.term[i].P[1,3]*system.branches.term[i].C[3] +
            system.branches.term[i].V0[3]);
            system.branches.term[i].P[1,4] = 5*mmHgToPa;
            # system.branches.term[i].P[1,4] = 4.9*mmHgToPa;
            system.branches.term[i].V[1,4] = (
            system.branches.term[i].P[1,4]*system.branches.term[i].C[4] +
            system.branches.term[i].V0[4]);
            # system.branches.term[i].P[1,5] = 4.5*mmHgToPa;
            system.branches.term[i].P[1,5] = 4.5*mmHgToPa;
            system.branches.term[i].V[1,5] = (
            system.branches.term[i].P[1,5]*system.branches.term[i].C[5] +
            system.branches.term[i].V0[5]);
            if system.branches.group[i] == "upper"
                numupper+=1
            end
        else
            # append!(system.branches.term.[i].P[:],NaN*[zeros(2,2)]);
            # append!(system.branches.term.[i].V[:],NaN*[zeros(2,2)]);
            # append!(system.branches.term.[i].Q[:],NaN*[zeros(2,2)]);
        end
    end

    if numupper > 0
        system.svc.P[1,1] = 4.5*mmHgToPa;
        system.svc.V[1,1] = system.svc.P[1,1]*system.svc.C + system.svc.V0;
        system.ivc.P[1,1] = 4.5*mmHgToPa;
        system.ivc.V[1,1] = system.ivc.P[1,1]*system.ivc.C + system.ivc.V0;
    else
        system.ivc.P[1,1] = 3*mmHgToPa;
        system.ivc.V[1,1] = system.ivc.P[1,1]*system.ivc.C + system.ivc.V0;
    end
end
