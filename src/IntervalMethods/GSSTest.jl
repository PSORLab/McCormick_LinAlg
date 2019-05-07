using IntervalArithmetic, LinearAlgebra
#Trying to decide how useful this search will be, timestep doesnt seem to impact det range over 0.
include("GoldenSectionSearch.jl")
include("GAUSSinv.jl")

Am = 100 * rand(10,10) .* rand([-1,1,1], 10,10)
r1 = .5
A_1 = map(x->Interval(x-r1, x+r1), Am)
I_ = fill(zero(Interval{Float64}), (10,10))
for i in 1:10
      I_[i,i] = Interval(1.0,1.0)
end
ts = [.0001, .001, .001, .01, .1, 1, 5, 10, 100]

detList = []

for t in ts
      push!(detList, gaussinvdet(I - A_1 * t))
end
detListZero = map(x->contains_zero(x), detList)
