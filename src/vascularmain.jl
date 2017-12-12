using CVModule
using MAT

filename = "arterytree.csv";

system = buildall(filename;numbeatstotal=1,restart="yes");

n = system.solverparams.nstart;

tic();
while system.solverparams.numbeats < system.solverparams.numbeatstotal
    predictorfluxes!(system,n);
    predictorstep!(system,n);
    correctorfluxes!(system,n);
    correctorstep!(system,n);
    applyendbcs!(system,n);
    splitinvariants!(system,n);
    solvesplits!(system,n);
    arterialpressure!(system,n);
    regulateall!(system,n);
    n+=1
end
toc()

updatevolumes!(system,n);

# file = matopen("solution.mat", "w")
# write(file, "system", system)
# close(file)
