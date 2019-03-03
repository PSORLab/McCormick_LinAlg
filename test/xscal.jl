#@test lines for XSCAL etc
#Meant to be the easiest test file to make
@testset "Test XSCAL" begin

#mctol = 1E-4
#= m = MC{2}(2.0, 3.0, IntervalType(1.0,4.0),
             seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false) =#
m = SVector{1,MC}(MC{3}(4.0, 5.0, IntervalType(4,5), SVector{3,Float64}([4.0, 5.0, 6.0]), SVector{3,Float64}([3.0, 2.0, 1.0]), false))

y = XSCAL(m, 2.0)[1]
yref = MC{3}(8.0, 10.0, IntervalType(8,10), SVector{3,Float64}([8.0, 10.0, 12.0]), SVector{3,Float64}([6.0, 4.0, 2.0]), false)

@test y.cv == yref.cv
@test y.cc == yref.cc
@test y.Intv.lo == yref.Intv.lo
@test y.Intv.hi == yref.Intv.hi
@test y.cv_grad[1] == yref.cv_grad[1]
@test y.cv_grad[2] == yref.cv_grad[2]
@test y.cv_grad[3] == yref.cv_grad[3]
@test y.cnst == yref.cnst

y = XSCAL(m, -2.0)[1]
yref = MC{3}(-10.0, -8.0, IntervalType(-10,-8), SVector{3,Float64}([-6.0, -4.0, -2.0]), SVector{3,Float64}([-8.0, -10.0, -12.0]), false)

@test y.cv == yref.cv
@test y.cc == yref.cc
@test y.Intv.lo == yref.Intv.lo
@test y.Intv.hi == yref.Intv.hi
@test y.cv_grad[1] == yref.cv_grad[1]
@test y.cv_grad[2] == yref.cv_grad[2]
@test y.cv_grad[3] == yref.cv_grad[3]
@test y.cnst == yref.cnst

end
