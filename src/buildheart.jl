type Activation
    a::Float64
    b::Float64
    k0::Float64
    k1::Float64
    th::Array{Float64,1}
    tce::Array{Float64,1}

    function Activation()
        this = new()
        this.a = 1.2;
        this.b = 0.3;
        this.k0 = 0.2;
        this.k1 = 0.2;
        this.th = [0.8];
        this.tce = [0.36];
        return this
    end
end

type LeftVentricle
    V0::Float64
    Emin::Float64
    V::Vector{Float64}
    P::Vector{Float64}
    E::Vector{Float64}
    Emax::Vector{Float64}

    function LeftVentricle()
        this = new()
        this.V0 = 10*cm3Tom3;
        this.Emin = 0.0283*mmHgToPa/cm3Tom3;
        this.Emax = [2*mmHgToPa/cm3Tom3];
        this.V = Vector{Float64}[];
        this.P = Vector{Float64}[];
        this.E = Vector{Float64}[];
        return this
    end
end

type LeftAtrium
    V0::Float64
    E::Float64
    R::Float64
    L::Float64
    V::Vector{Float64}
    P::Vector{Float64}
    Q::Vector{Float64}


    function LeftAtrium()
        this = new()
        this.V0 = 10*cm3Tom3;
        this.R = 3.6e-3*mmHgToPa/cm3Tom3;
        this.L = 3e-5*mmHgToPa/cm3Tom3;
        this.E = 0.13*mmHgToPa/cm3Tom3;
        this.V = Vector{Float64}[];
        this.P = Vector{Float64}[];
        this.Q = Vector{Float64}[];
        return this
    end
end

type RightVentricle
    V0::Float64
    Emin::Float64
    L::Float64
    V::Vector{Float64}
    P::Vector{Float64}
    E::Vector{Float64}
    Q::Vector{Float64}
    Emax::Vector{Float64}

    function RightVentricle()
        this = new()
        this.V0 = 10*cm3Tom3;
        this.Emin = 0.0283*mmHgToPa/cm3Tom3;
        this.Emax = [0.36*mmHgToPa/cm3Tom3];
        this.L = 2.16e-4*mmHgToPa/cm3Tom3;
        this.V = Vector{Float64}[];
        this.P = Vector{Float64}[];
        this.E = Vector{Float64}[];
        this.Q = Vector{Float64}[];
        return this
    end
end

type RightAtrium
    V0::Float64
    E::Float64
    R::Float64
    L::Float64
    V::Vector{Float64}
    P::Vector{Float64}
    Q::Vector{Float64}


    function RightAtrium()
        this = new()
        this.V0 = 10*cm3Tom3;
        this.R = 4.85e-3*mmHgToPa/cm3Tom3;
        this.L = 5e-5*mmHgToPa/cm3Tom3;
        this.E = 0.16*mmHgToPa/cm3Tom3;
        this.V = Vector{Float64}[];
        this.P = Vector{Float64}[];
        this.Q = Vector{Float64}[];
        return this
    end
end

type AorticValve
    Kvo::Float64
    Kvc::Float64
    leff::Float64
    Po::Float64
    Pc::Float64
    Aann::Float64
    Ks::Float64
    zeta::Vector{Float64}

    function AorticValve()
        this = new()
        this.Kvo = 0.1;
        this.Kvc = 0.1;
        this.leff = 0.01;
        this.Po = 0;
        this.Pc = 0;
        this.Aann = 6.8e-4;
        this.Ks = 4e-9/cm3Tom3;
        this.zeta = Vector{Float64}[];
        return this
    end
end

type Heart
    activation::Activation
    lv::LeftVentricle
    la::LeftAtrium
    rv::RightVentricle
    ra::RightAtrium
    av::AorticValve

    function Heart()
        this = new()
        this.activation = Activation();
        this.lv = LeftVentricle();
        this.la = LeftAtrium();
        this.rv = RightVentricle();
        this.ra = RightAtrium();
        this.av = AorticValve();
        return this
    end
end
