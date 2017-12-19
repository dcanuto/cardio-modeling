type CNS # central nervous system
    # sympathetic, parasympathetic activations
    np::Vector{Float64}
    ns::Vector{Float64}
    # heart rate and maximum ventriclular elastance
    H::Vector{Float64}
    Emaxlv::Vector{Float64}
    Emaxrv::Vector{Float64}
    # lower/upper body compartmental resistances, compliances, unstressed vols.
    R2L::Vector{Float64}
    R3L::Vector{Float64}
    R2U::Vector{Float64}
    R3U::Vector{Float64}
    C4L::Vector{Float64}
    C5L::Vector{Float64}
    C4U::Vector{Float64}
    C5U::Vector{Float64}
    V4L::Vector{Float64}
    V5L::Vector{Float64}
    V4U::Vector{Float64}
    V5U::Vector{Float64}
    # feedback model weights (alpha:ns, beta:np, gamma:bias)
    alphaH::Float64
    betaH::Float64
    gammaH::Float64
    alphalv::Float64
    gammalv::Float64
    alpharv::Float64
    gammarv::Float64
    alphaR2L::Float64
    gammaR2L::Float64
    alphaR3L::Float64
    gammaR3L::Float64
    alphaR2U::Float64
    gammaR2U::Float64
    alphaR3U::Float64
    gammaR3U::Float64
    alphaC4L::Float64
    gammaC4L::Float64
    alphaC5L::Float64
    gammaC5L::Float64
    alphaC4U::Float64
    gammaC4U::Float64
    alphaC5U::Float64
    gammaC5U::Float64
    alphaV4L::Float64
    gammaV4L::Float64
    alphaV5L::Float64
    gammaV5L::Float64
    alphaV4U::Float64
    gammaV4U::Float64
    alphaV5U::Float64
    gammaV5U::Float64
    # time constants
    tauH::Float64
    tauEmax::Float64
    tauC::Float64
    tauR::Float64
    tauV::Float64
    # activation function parameters
    Ptarget::Float64
    Paverage::Vector{Float64}
    baselinefirerate::Float64
    a::Float64
    b::Float64
    c::Float64
    # indices for baroreflex pressure averaging
    avgindexstart::Int32
    avgindexend::Int32


    function CNS()
        this = new()
        this.ns = Vector{Float64}[];
        this.np = Vector{Float64}[];
        this.H = Vector{Float64}[];
        this.Emaxlv = Vector{Float64}[];
        this.Emaxrv = Vector{Float64}[];
        this.R2L = Vector{Float64}[];
        this.R3L = Vector{Float64}[];
        this.R2U = Vector{Float64}[];
        this.R3U = Vector{Float64}[];
        this.C4L = Vector{Float64}[];
        this.C5L = Vector{Float64}[];
        this.C4U = Vector{Float64}[];
        this.C5U = Vector{Float64}[];
        this.V4L = Vector{Float64}[];
        this.V5L = Vector{Float64}[];
        this.V4U = Vector{Float64}[];
        this.V5U = Vector{Float64}[];
        this.tauH = 2;
        this.tauEmax = 2;
        this.tauC = 20;
        this.tauR = 6;
        this.tauV = 20;
        this.Ptarget = 96.5*mmHgToPa;
        this.Paverage = Vector{Float64}[];
        this.avgindexstart = 1;
        this.avgindexend = 0;
        this.baselinefirerate = 0.25;
        this.a = 1;
        this.b = -log(this.baselinefirerate);
        this.c = 6;
        return this
    end
end
