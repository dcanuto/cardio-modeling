module CVModule # All functions/types needed for closed-loop CV system solver

using NumericalIntegration
using MAT
using Interpolations

# conversion factors and solver parameters
include("conversions.jl")
include("solverparams.jl")

# function for loading 1D domain empirical data from CSV file
include("loadtexttree.jl")

# type definitions
include("buildbranches.jl")
include("buildheart.jl")
include("buildvenacava.jl")
include("buildlungs.jl")
include("buildcns.jl")
include("buildhemo.jl")
include("buildall.jl")

# memory allocators for solution variables
include("calcbranchprops.jl")
include("discretizebranches.jl")
include("assignterminals.jl")
include("discretizeperiphery.jl")
include("discretizeheart.jl")
include("discretizelungs.jl")
include("discretizecns.jl")
include("updatediscretization.jl")

# initial conditions for each submodel
include("applybranchics.jl")
include("applyperipheryics.jl")
include("applyheartics.jl")
include("applylungics.jl")
include("applycnsics.jl")
include("applycustomics.jl")

# heart model
include("elastancefn.jl")
include("elastancemodel.jl")
include("setnumbeats.jl")

# volume tracker, error corrector for mass conservation
include("updatevolumes.jl")
include("correctvolume.jl")

# Lax-Wendroff method for interior points of arterial network
include("predictorfluxes.jl")
include("predictorstep.jl")
include("correctorfluxes.jl")
include("correctorstep.jl")

# 0D-1D coupling and 0D updates
include("applyendbcs.jl")
include("coupledistal.jl")
include("endinvariants.jl")
include("updateterms.jl")
include("updatevc.jl")
include("updaterh.jl")
include("updatelungs.jl")
include("updatela.jl")
include("coupleproximal.jl")
include("rootinvariant.jl")
include("fav.jl")
include("Jav.jl")
include("newtonav.jl")
include("linesearch.jl")
include("fdist.jl")
include("Jdist.jl")
include("newtondist.jl")
include("linedist.jl")

# 1D interior junction updates
include("splitinvariants.jl")
include("solvesplits.jl")
include("newton.jl")
include("fsingle.jl")
include("fdouble.jl")
include("ftriple.jl")
include("fquad.jl")
include("Jsingle.jl")
include("Jdouble.jl")
include("Jtriple.jl")
include("Jquad.jl")

# 1D pressure updates
include("arterialpressure.jl")

# regulation
include("regulateall.jl")
include("regulateperiphery.jl")
include("reflexpressure.jl")
include("cnsactivations.jl")
include("regulateheart.jl")

# hemorrhage/tourniquet application
include("applytourniquet.jl")
include("applyhemoics.jl")
include("modelhemo.jl")

export loadtexttree
export calcbranchprops!
export discretizebranches!
export assignterminals!
export discretizeperiphery!
export discretizeheart!
export discretizelungs!
export discretizecns!
export applybranchics!
export applyperipheryics!
export applyheartics!
export applylungics!
export applycnsics!
export applycustomics!
export elastancefn!
export elastancemodel!
export updatevolumes!
export buildall
export predictorfluxes!
export predictorstep!
export correctorfluxes!
export correctorstep!
export applyendbcs!
export coupledistal!
export endinvariants!
export updateterms!
export updatevc!
export updaterh!
export updatelungs!
export updatela!
export coupleproximal!
export rootinvariant!
export splitinvariants!
export solvesplits!
export newton!
export fsingle
export fdouble
export ftriple
export fquad
export Jsingle
export Jdouble
export Jtriple
export Jquad
export arterialpressure!
export regulateall!
export regulateperiphery!
export setnumbeats!
export reflexpressure!
export cnsactivations!
export regulateheart!
export updatediscretization!
export applytourniquet!

end
