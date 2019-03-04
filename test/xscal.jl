#@test lines for XSCAL etc
#Meant to be the easiest test file to make
@testset "Test XSCAL" begin

#mctol = 1E-4
#= m = MC{2}(2.0, 3.0, IntervalType(1.0,4.0),
             seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false) =#
m = SVector{2,MC}(MC{3}(4.0, 5.0, IntervalType(4,5), SVector{3,Float64}([4.0, 5.0, 6.0]), SVector{3,Float64}([3.0, 2.0, 1.0]), false),MC{3}(4.0, 5.0, IntervalType(4,5), SVector{3,Float64}([4.0, 5.0, 6.0]), SVector{3,Float64}([3.0, 2.0, 1.0]), false))

y = XSCAL(m, 2.0)
y1 = y[1]
y2 = y[2]
yref = MC{3}(8.0, 10.0, IntervalType(8,10), SVector{3,Float64}([8.0, 10.0, 12.0]), SVector{3,Float64}([6.0, 4.0, 2.0]), false)

@test y1.cv == yref.cv
@test y1.cc == yref.cc
@test y1.Intv.lo == yref.Intv.lo
@test y1.Intv.hi == yref.Intv.hi
@test y1.cv_grad[1] == yref.cv_grad[1]
@test y1.cv_grad[2] == yref.cv_grad[2]
@test y1.cv_grad[3] == yref.cv_grad[3]
@test y1.cnst == yref.cnst

@test y2.cv == yref.cv
@test y2.cc == yref.cc
@test y2.Intv.lo == yref.Intv.lo
@test y2.Intv.hi == yref.Intv.hi
@test y2.cv_grad[1] == yref.cv_grad[1]
@test y2.cv_grad[2] == yref.cv_grad[2]
@test y2.cv_grad[3] == yref.cv_grad[3]
@test y2.cnst == yref.cnst

y = XSCAL(m, -2.0)
y1 = y[1]
y2 = y[2]
yref = MC{3}(-10.0, -8.0, IntervalType(-10,-8), SVector{3,Float64}([-6.0, -4.0, -2.0]), SVector{3,Float64}([-8.0, -10.0, -12.0]), false)

@test y1.cv == yref.cv
@test y1.cc == yref.cc
@test y1.Intv.lo == yref.Intv.lo
@test y1.Intv.hi == yref.Intv.hi
@test y1.cv_grad[1] == yref.cv_grad[1]
@test y1.cv_grad[2] == yref.cv_grad[2]
@test y1.cv_grad[3] == yref.cv_grad[3]
@test y1.cnst == yref.cnst

@test y2.cv == yref.cv
@test y2.cc == yref.cc
@test y2.Intv.lo == yref.Intv.lo
@test y2.Intv.hi == yref.Intv.hi
@test y2.cv_grad[1] == yref.cv_grad[1]
@test y2.cv_grad[2] == yref.cv_grad[2]
@test y2.cv_grad[3] == yref.cv_grad[3]
@test y2.cnst == yref.cnst

end
