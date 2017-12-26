type Hemorrhage
    hID::Array{Int64,1}
    tID::Array{Int64,1}
    Amax::Array{Float64,1}
    Amin::Array{Float64,1}
    bmax::Array{Float64,1}
    bmin::Array{Float64,1}
    Vloss::Float64
    Vlossinit::Float64
    totalloss::Float64
    Ph::Float64
    ttn::Float64
    injured::Vector{Bool}

    function Hemorrhage()
        this = new()
        this.hID = [41];
        this.tID = [39,40];
        this.Amax = Array{Float64,1}[];
        this.Amin = Array{Float64,1}[];
        this.bmax = Array{Float64,1}[];
        this.bmin = Array{Float64,1}[];
        this.Vloss = 0;
        this.Vlossinit = 0;
        this.totalloss = 550*cm3Tom3;
        this.Ph = 0*mmHgToPa;
        this.ttn = 10;
        this.injured = Vector{String}[];
        return this
    end
end
