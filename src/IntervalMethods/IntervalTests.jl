#Yes this isnt in the test file. Just while this folder is in development
using IntervalArithmetic
include("SpecialFunctions.jl")
include("CRAM.jl")
include("CRAMinv.jl")
include("GAUSS.jl")
include("GAUSSinv.jl")
include("GERSCH.jl")
include("GERSCHinv.jl")
Am = [1. 2. 3.
      4. 6. 7.
      5. 9. 8.] #Midpoint/Center Matrix
r1 = .1
r2 = .01
A_1 = map(x->Interval(x-r1, x+r1), Am) #Interval Matrix
A_2 = map(x->Interval(x-r2, x+r2), Am) #Interval Matrix
HULL = [Interval(4.060, 14.880), Interval(8.465, 9.545)]#From paper
DET = det(Am) #Actual det of Am is 9
data = Array{Any,2}(undef,(8,3))
data[1,:] = ["method", "r=0.1", "r=0.01"]
data[2:end, 1] = ["HULL";"GE"; "GEinv"; "CRAM"; "CRAMinv"; "GERSCH"; "GERSCHinv"]

data[2:end, 2] = [HULL[1]; gaussdet(A_1)[1]; gaussinvdet(A_1)[1]; cramdet(A_1); craminvdet(A_1); GERSCH(A_1); GERSCHinv(A_1)]
data[2:end, 3] = [HULL[2]; gaussdet(A_2)[1]; gaussinvdet(A_2)[1]; cramdet(A_2); craminvdet(A_2); GERSCH(A_2); GERSCHinv(A_2)]
display(data)
#Add GERSCH, HADinv, and make DataFrame?
#CRAM doesnt seem to be working, bounds do not contain 9.
#GAUSS and GAUSSinv still have to methods in them to compare
display(["method" "r = 0.1" "r = 0.01"
"HULL" "[4.060, 14.880]" "[8.465, 9.545]"
"GE" "[3.000, 21.857]" "[8.275, 9.789]"
"GEinv" "[3.600, 18.000]" "[8.460, 9.560]"
"GElu" "[1.440, 22.482]" "[8.244, 9.791]"
"CRAM" "[-∞, ∞]" "[8.326, 9.765]"
"CRAMinv" "[3.594, 78.230]" "[8.460, 9.588]"
"CRAMlu" "[-∞, ∞]" "[8.244, 9.863]"
"HAD" "[-526.712, 526.712]" "[-493.855, 493.855]"
"HADinv" "[-16.801, 16.801]" "[-9.563, 9.563]"])
