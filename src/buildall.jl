type CVSystem # entire solution
    heart::Heart
    branches::ArterialBranches
    svc::VenaCava
    ivc::VenaCava
    lungs::Lungs
    cns::CNS
    solverparams::SolverParams
    t::Vector{Float64}
    arterialvolume::Float64
    peripheralvolume::Float64
    vcvolume::Float64
    heartvolume::Float64
    lungvolume::Float64
    initialvolume::Float64
    finalvolume::Float64

    function CVSystem(filename="test.csv")
        this = new()
        this.heart = Heart();
        this.branches = ArterialBranches(filename);
        this.svc = VenaCava();
        this.ivc = VenaCava();
        this.lungs = Lungs();
        this.cns = CNS();
        this.solverparams = SolverParams();
        this.t = Vector{Float64}[];
        return this
    end
end

# build solution struct
function buildall(filename="test.csv";numbeatstotal=0,restart="yes")
    system = CVSystem(filename);
    system.solverparams.numbeatstotal = numbeatstotal;
    calcbranchprops!(system)
    discretizebranches!(system);
    assignterminals!(system);
    discretizeperiphery!(system);
    discretizeheart!(system);
    discretizelungs!(system);
    discretizecns!(system);
    applybranchics!(system);
    applyperipheryics!(system);
    applyheartics!(system);
    applylungics!(system);
    applycnsics!(system);
    applycustomics!(system);
    updatevolumes!(system,0);
    return system
end
