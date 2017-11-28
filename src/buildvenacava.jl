type VenaCava
    R::Float64
    L::Float64
    C::Float64
    V0::Float64
    P::Vector{Float64}
    V::Vector{Float64}
    Q::Vector{Float64}

    function VenaCava()
        this = new()
        this.P = Vector{Float64}[];
        this.V = Vector{Float64}[];
        this.Q = Vector{Float64}[];
        return this
    end
end
