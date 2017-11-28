type ArterialTerminal # lumped peripheral vasculature
    C::Vector{Float64}
    R::Vector{Float64}
    L::Vector{Float64}
    V0::Vector{Float64}
    P::Any
    V::Any
    Q::Any

    function ArterialTerminal()
        this = new()
        this.C = Vector{Float64}[];
        this.R = Vector{Float64}[];
        this.L = Vector{Float64}[];
        this.V0 = Vector{Float64}[];
        this.P = Array{Float64,2}[];
        this.V = Array{Float64,2}[];
        this.Q = Array{Float64,2}[];
        return this
    end
end

type ArterialBranches # 1D arterial domain
    name::Vector{String}
    parentname::Vector{String}
    ID::Vector{Int64}
    parentID::Vector{Int64}
    lengthincm::Vector{Float64}
    radiusincm::Vector{Float64}
    thicknessincm::Vector{Float64}
    YoungsModinMPa::Vector{Float64}
    A0::Any
    beta::Any
    c0::Any
    k::Vector{Float64}
    group::Vector{String}
    A::Any
    Q::Any
    P::Any
    Fp::Any
    Fbarforward::Any
    Fbarbackward::Any
    Abackward::Any
    Aforward::Any
    Qbackward::Any
    Qforward::Any
    W1end::Vector{Float64}
    W1::Vector{Float64}
    W2::Vector{Float64}
    W1root::Float64
    W2root::Float64
    children::Any
    term::Any
    Rao::Float64

    function ArterialBranches(filename="test.csv")
        this = new()
        this.name = Vector{String}[];
        this.parentname = Vector{String}[];
        this.ID = Vector{Int64}[];
        this.parentID = Vector{Int64}[];
        this.children = Vector{Int64}[];
        this.lengthincm = Vector{Float64}[];
        this.radiusincm = Vector{Float64}[];
        this.thicknessincm = Vector{Float64}[];
        this.YoungsModinMPa = Vector{Float64}[];
        this.A0 = Array{Float64,1}[];
        this.beta = Array{Float64,1}[];
        this.c0 = Array{Float64,1}[];
        this.k = Vector{Float64}[];
        this.group = Vector{String}[];
        this.A = Array{Float64,2}[];
        this.Q = Array{Float64,2}[];
        this.P = Array{Float64,2}[];
        this.Fp = Array{Float64,1}[];
        this.Fbarforward = Array{Float64,1}[];
        this.Fbarbackward = Array{Float64,1}[];
        this.Abackward = Array{Float64,1}[];
        this.Aforward = Array{Float64,1}[];
        this.Qbackward = Array{Float64,1}[];
        this.Qforward = Array{Float64,1}[];
        this.W1end = Vector{Float64}[];
        this.W1 = Vector{Float64}[];
        this.W2 = Vector{Float64}[];
        this.W1root = Float64(0);
        this.W2root = Float64(0);
        this.Rao = 0.01*mmHgToPa/cm3Tom3;

        # pull in artery data from text file
        temp = loadtexttree(filename);

        foreach((x)->push!(this.name,get(x)),temp[:Name])
        foreach((x)->push!(this.parentname,get(x)),temp[:ParentName])
        foreach((x)->push!(this.ID,get(x)),temp[:ID])
        foreach((x)->push!(this.parentID,get(x)),temp[:parentID])
        for i in 1:size(temp[:children_1],1)
            push!(this.group,get(temp[i,:group],"NA"))
            push!(this.children,[get(temp[i,:children_1],0)])
            push!(this.children[i],get(temp[i,:children_2],0))
            push!(this.children[i],get(temp[i,:children_3],0))
            push!(this.children[i],get(temp[i,:children_4],0))
            deleteat!(this.children[i],findin(this.children[i],0))
        end
        foreach((x)->push!(this.lengthincm,get(x)),temp[:Length_cm])
        foreach((x)->push!(this.radiusincm,get(x)),temp[:Radius_cm])
        foreach((x)->push!(this.thicknessincm,get(x)),temp[:Thickness_cm])
        foreach((x)->push!(this.YoungsModinMPa,get(x)),temp[:YoungsModulus_MPa])

        this.term = Vector{ArterialTerminal}(length(this.ID));

        return this
    end
end
