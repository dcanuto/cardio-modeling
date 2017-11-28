# cardio-modeling
Closed-loop model of the cardiovascular system with autonomic regulation. No fancy frills for interactive sessions yet (i.e., running a case requires manual changes to the `vascular_main.jl` script). First, open a Julia session in the `src\` directory, then build the necessary types and functions:

    include("CVModule.jl")  

To run a case for `n` cardiac cycles, open the `vascular_main.jl` script and alter line 5:

    system = buildall(filename;numbeatstotal=n,restart="yes");
