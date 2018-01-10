using CVModule
using MAT

function vascularmain()

# filename = "test.csv";
filename = "uc11.mat";
rstflag = "yes"
hemoflag = "yes"
saveflag = "yes"

system = buildall(filename;numbeatstotal=1,restart=rstflag,injury=hemoflag);

n = system.solverparams.nstart;

tic();
while system.solverparams.numbeats < system.solverparams.numbeatstotal
    predictorfluxes!(system,n);
    predictorstep!(system,n);
    correctorfluxes!(system,n);
    correctorstep!(system,n);
    applyendbcs!(system,n);
    splitinvariants!(system,n);
    if hemoflag == "no"
        solvesplits!(system,n);
    elseif hemoflag == "yes"
        solvesplits!(system,n,hemoflag);
        applytourniquet!(system,n);
    end
    arterialpressure!(system,n);
    regulateall!(system,n);
    n+=1
end
toc()

updatevolumes!(system,n);

if saveflag == "yes"
    file = matopen("uc12.mat", "w")
    write(file, "system", system)
    close(file)
end

return system, n

end
