type Activation
    th::Array{Float64,1}
    tau1::Float64
    tau2::Float64
    m1::Float64
    m2::Float64
    k::Array{Float64,1}

    function Activation(old=Dict("a"=>0),restart="no")
        this = new()
        this.m1 = 1.32;
        this.m2 = 27.4;
        if restart == "no"
            this.th = [0.8];
            this.tau1 = 0.269*this.th[1];
            this.tau2 = 0.452*this.th[1];
            t = linspace(0,this.th[1],10000);
            g1 = (t/this.tau1).^this.m1;
            g2 = (t/this.tau2).^this.m2;
            this.k = [maximum((g1./(1+g1)).*(1./(1+g2)))^-1];
        elseif restart == "yes"
            this.th = [(old["th"])[end]];
            this.tau1 = old["tau1"];
            this.tau2 = old["tau2"];
            this.k = [old["k"][end]];
        end
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

    function LeftVentricle(old=Dict("a"=>0),restart="no")
        this = new()
        this.V0 = 10*cm3Tom3;
        this.Emin = 0.0283*mmHgToPa/cm3Tom3;
        if restart == "no"
            this.Emax = [3*mmHgToPa/cm3Tom3];
        elseif restart == "yes"
            this.Emax = [old["Emax"][end]];
        end
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

    function RightVentricle(old=Dict("a"=>0),restart="no")
        this = new()
        this.V0 = 10*cm3Tom3;
        this.Emin = 0.0283*mmHgToPa/cm3Tom3;
        if restart == "no"
            this.Emax = [0.4*mmHgToPa/cm3Tom3];
        elseif restart == "yes"
            this.Emax = [old["Emax"][end]];
        end
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
        this.Kvo = 0.16;
        this.Kvc = 0.16;
        this.leff = 0.01;
        this.Po = 0;
        this.Pc = 0;
        this.Aann = 3.8e-4;
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

    function Heart(old=Dict("a"=>0),restart="no")
        this = new()
        if restart == "yes"
            act = old["activation"];
            lv = old["lv"];
            rv = old["rv"];
            this.activation = Activation(act,restart);
            this.lv = LeftVentricle(lv,restart);
            this.rv = RightVentricle(rv,restart);
        elseif restart == "no"
            this.activation = Activation();
            this.lv = LeftVentricle();
            this.rv = RightVentricle();
        end
        this.la = LeftAtrium();
        this.ra = RightAtrium();
        this.av = AorticValve();
        return this
    end
end
