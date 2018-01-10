type CVSystem # entire solution
    heart::Heart
    branches::ArterialBranches
    svc::VenaCava
    ivc::VenaCava
    lungs::Lungs
    hemo::Hemorrhage
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

    function CVSystem(filename="test.csv",restart="no")
        this = new()
        if restart == "no"
            this.heart = Heart();
            this.branches = ArterialBranches(filename);
        elseif restart == "yes"
            vars = matread(filename);
            sys = vars["system"];
            heart = sys["heart"];
            branches = sys["branches"];
            this.heart = Heart(heart,restart)
            this.branches = ArterialBranches(filename,branches,restart)
        else
            error("Keyword restart must either be yes or no. Aborting.")
        end
        this.svc = VenaCava();
        this.ivc = VenaCava();
        this.lungs = Lungs();
        this.cns = CNS();
        this.hemo = Hemorrhage();
        this.solverparams = SolverParams();
        this.t = Vector{Float64}[];
        return this
    end
end

# build solution struct
function buildall(filename="test.csv";numbeatstotal=1,restart="no",injury="no")
    if restart == "no"
        system = CVSystem(filename);
        system.solverparams.numbeatstotal = numbeatstotal;
        calcbranchprops!(system);
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
    elseif restart == "yes"
        vars = matread(filename);
        sys = vars["system"];
        branches = sys["branches"];
        term = branches["term"];
        heart = sys["heart"];
        lungs = sys["lungs"];
        cns = sys["cns"];
        system = CVSystem(filename,restart);
        system.solverparams.numbeatstotal = numbeatstotal;
        calcbranchprops!(system,branches,restart);
        discretizebranches!(system,sys,restart);
        assignterminals!(system,term,restart);
        discretizeperiphery!(system);
        discretizeheart!(system);
        discretizelungs!(system);
        discretizecns!(system);
        applybranchics!(system,sys,restart);
        applyperipheryics!(system,sys,restart);
        applyheartics!(system,heart,restart);
        applylungics!(system,lungs,restart);
        applycnsics!(system,cns,restart);
        applyhemoics!(system,sys);
    end
    updatevolumes!(system,0);
    return system
end
