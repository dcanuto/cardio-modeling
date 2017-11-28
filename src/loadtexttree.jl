using CSV

function loadtexttree(filename::String)
    branches = CSV.read(filename,null="");
    return branches
end
