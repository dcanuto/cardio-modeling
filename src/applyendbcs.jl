function applyendbcs!(system::CVSystem,n::Int64)
    # update periphery
    coupledistal!(system,n);
    updateterms!(system,n);
    updatevc!(system,n);

    # update right heart
    updaterh!(system,n);

    # update pulmonary circulation
    updatelungs!(system,n);

    # update left heart
    updatela!(system,n);
    coupleproximal!(system,n);

    # update heart chamber pressures
    elastancemodel!(system,n+1);
end
