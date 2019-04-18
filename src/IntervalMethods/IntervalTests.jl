#Yes this isnt in the test file. Just while this folder is in development
using IntervalArithmetic
include("CRAMinv.jl")
include("GaussElimDeterminant.jl")
include("GaussElimInvDeterminant.jl")
Am = [1. 2. 3.
      4. 6. 7.
      5. 9. 8.] #Midpoint/Center Matrix
r1 = .1
r2 = .01
A_1 = map(x->Interval(x-r1, x+r1), Am) #Interval Matrix
A_2 = map(x->Interval(x-r2, x+r2), Am) #Interval Matrix
HULL = [[4.060, 14.880] [8.465, 9.545]]#From paper
DET = det(Am) #Actual det of Am is 9
data = Array{Any,2}(undef,(4,3))
data[1,:] = ["method", "r=0.1", "r=0.01"]
data[2:end, 1] = ["GE"; "GEinv"; "CRAMinv"]

data[2:end, 2] = [gaussdet(A_1); gaussinvdet(A_1); CRAMinvdet(A_1) ]
data[2:end, 3] = [gaussdet(A_2); gaussinvdet(A_2); CRAMinvdet(A_2) ]
