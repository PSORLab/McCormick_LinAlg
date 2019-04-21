 #replica of table 1 from DetofIntMatrices Paper (HORA´CEK)
using IntervalArithmetic, LinearAlgebra
include("SpecialFunctions.jl")
include("CRAM.jl")
include("CRAMinv.jl")
include("GAUSS.jl")
include("GAUSSinv.jl")
include("GERSCH.jl")
include("GERSCHinv.jl")
include("HAD.jl")
include("HADinv.jl")
Am = [1. 2. 3.
      4. 6. 7.
      5. 9. 8.] #Midpoint/Center Matrix
r1 = .1
r2 = .01
A_1 = map(x->Interval(x-r1, x+r1), Am) #Interval Matrix with radius =.1
A_2 = map(x->Interval(x-r2, x+r2), Am) #Interval Matrix with radius =.01
HULL = [Interval(4.060, 14.880), Interval(8.465, 9.545)]#From paper
DET = det(Am) #det of Am is 9

data = Array{Any,2}(undef,(9,3))
data[1,:] = ["method", "r=0.1", "r=0.01"]
data[2:end, 1] = ["HULL";"GE"; "GEinv"; "CRAM"; "CRAMinv"; "GERSCH"; "GERSCHinv"; "HADinv"]

data[2:end, 2] = [HULL[1]; gaussdet(A_1)[1]; gaussinvdet(A_1)[1]; cramdet(A_1); craminvdet(A_1); GERSCH(A_1); GERSCHinv(A_1); HADinv(A_1)]
data[2:end, 3] = [HULL[2]; gaussdet(A_2)[1]; gaussinvdet(A_2)[1]; cramdet(A_2); craminvdet(A_2); GERSCH(A_2); GERSCHinv(A_2); HADinv(A_2)]
println("Our results")
display(data);
#GAUSS and GAUSSinv still have to methods in them to compare
println("results from (HORA´CEK)")
display(["method" "r = 0.1" "r = 0.01"
"HULL" "[4.060, 14.880]" "[8.465, 9.545]"
"GE" "[3.000, 21.857]" "[8.275, 9.789]"
"GEinv" "[3.600, 18.000]" "[8.460, 9.560]"
"GElu" "[1.440, 22.482]" "[8.244, 9.791]"
"CRAM" "[-∞, ∞]" "[8.326, 9.765]"
"CRAMinv" "[3.594, 78.230]" "[8.460, 9.588]"
"CRAMlu" "[-∞, ∞]" "[8.244, 9.863]"
"HAD" "[-526.712, 526.712]" "[-493.855, 493.855]"
"HADinv" "[-16.801, 16.801]" "[-9.563, 9.563]"
"GERSCH" "[-3132.927, 11089.567]" "[-2926.485, 10691.619]"
"GERSCHinv" "[-0.000, 72.000]" "[6.561, 11.979]"  ]);
