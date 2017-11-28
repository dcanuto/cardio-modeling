type Lungs
    Rp::Float64
    Ra::Vector{Float64}
    Rv::Vector{Float64}
    Ca::Vector{Float64}
    Cv::Vector{Float64}
    La::Float64
    Lv::Float64
    V0::Vector{Float64}
    Pp::Vector{Float64}
    Pa::Any
    Va::Any
    Qa::Any
    Pv::Any
    Vv::Any
    Qv::Any

    function Lungs()
        this = new()
        this.Rp = 0.0251*mmHgToPa/cm3Tom3
        this.Ra = [0.0227*mmHgToPa/cm3Tom3];
        push!(this.Ra,0.03*mmHgToPa/cm3Tom3);
        push!(this.Ra,0.021*mmHgToPa/cm3Tom3);
        this.Rv = [0.010*mmHgToPa/cm3Tom3];
        push!(this.Rv,0.010*mmHgToPa/cm3Tom3);
        this.Ca = [2.222*cm3Tom3/mmHgToPa];
        push!(this.Ca,1.481*cm3Tom3/mmHgToPa);
        push!(this.Ca,1.778*cm3Tom3/mmHgToPa);
        this.Cv = [13*cm3Tom3/mmHgToPa];
        push!(this.Cv,74*cm3Tom3/mmHgToPa);
        this.La = 5e-5*mmHgToPa/cm3Tom3;
        this.Lv = 5e-5*mmHgToPa/cm3Tom3;
        this.V0 = [50*cm3Tom3];
        append!(this.V0,[30,53,75,75]*cm3Tom3);
        this.Pp = Vector{Float64}[];
        this.Pa = Array{Float64,2}[];
        this.Va = Array{Float64,2}[];
        this.Qa = Array{Float64,2}[];
        this.Pv = Array{Float64,2}[];
        this.Vv = Array{Float64,2}[];
        this.Qv = Array{Float64,2}[];
        return this
    end
end
