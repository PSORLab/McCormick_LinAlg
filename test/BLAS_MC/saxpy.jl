#@test lines for XSCAL etc
#Meant to be the easiest test file to make
@testset "Test SAXPY" begin
#RELAXATIONS INSIDE THE INTERVAL BOUNDS
mctol = 1E-4
#= m = MC{2}(2.0, 3.0, IntervalType(1.0,4.0),
             seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false) =#
#To do make test vectors less similar
m1= MC{3}(4.0, 5.0, IntervalType(4.,5.), SVector{3,Float64}(4.0, 5.0, 6.0), SVector{3,Float64}(3.0, 2.0, 1.0), false)
m2= MC{3}(4.0, 5.0, IntervalType(4.,5.), SVector{3,Float64}(4.0, 5.0, 6.0), SVector{3,Float64}(3.0, 2.0, 1.0), false)
m3= MC{3}(4.0, 5.0, IntervalType(4., 5.), SVector{3,Float64}(4.0, 5.0, 6.0), SVector{3,Float64}(3.0, 2.0, 1.0),false)
m4= MC{3}(3.0, 4.0, IntervalType(3., 4.), SVector{3,Float64}(4.0, 5.0, 6.0), SVector{3,Float64}(3.0, 2.0, 1.0), false)
X = [m1,m2]
Y = [m3, m4]

a::Float64 = 6.3

y = SAXPY(a, X, Y)
y1 = y[1]
y2 = y[2]

yref = [MC{3}(29.2, 36.5, IntervalType(29.1999, 36.5), SVector{3,Float64}(29.2, 36.5, 43.8), SVector{3,Float64}(21.9, 14.6, 7.3), false),
                      MC{3}(28.2, 35.5, IntervalType(28.1999, 35.5), SVector{3,Float64}(29.2, 36.5, 43.8), SVector{3,Float64}(21.9, 14.6, 7.3), false)]
yref1 = yref[1]
yref2 = yref[2]

@test isapprox(y1.cv,      yref1.cv, atol=mctol)#lots of rounding error 21.9999 instead of 22.0
@test isapprox(y1.cc,      yref1.cc, atol=mctol)
@test isapprox(y1.Intv.lo, yref1.Intv.lo, atol=mctol)
@test isapprox(y1.Intv.hi, yref1.Intv.hi, atol=mctol)
@test isapprox(y1.cv_grad[1], yref1.cv_grad[1], atol=mctol)
@test isapprox(y1.cv_grad[2], yref1.cv_grad[2], atol=mctol)
@test isapprox(y1.cv_grad[3], yref1.cv_grad[3], atol=mctol)
@test isapprox(y1.cnst,    yref1.cnst, atol=mctol)

@test isapprox(y2.cv,      yref2.cv, atol=mctol)
@test isapprox(y2.cc,      yref2.cc, atol=mctol)
@test isapprox(y2.Intv.lo, yref2.Intv.lo, atol=mctol)
@test isapprox(y2.Intv.hi, yref2.Intv.hi, atol=mctol)
@test isapprox(y2.cv_grad[1], yref2.cv_grad[1], atol=mctol)
@test isapprox(y2.cv_grad[2], yref2.cv_grad[2], atol=mctol)
@test isapprox(y2.cv_grad[3], yref2.cv_grad[3], atol=mctol)
@test isapprox(y2.cnst,    yref2.cnst, atol=mctol)

a = -6.3
y = SAXPY(a, X, Y)
y1 = y[1]
y2 = y[2]

yref =[MC{3}(-27.5, -20.2, IntervalType(-27.5, -20.1999), SVector{3,Float64}(-14.9, -7.6, -0.3), SVector{3,Float64}(-22.2, -29.5, -36.8), false),
                    MC{3}(-28.5, -21.2, IntervalType(-28.5, -21.1999), SVector{3,Float64}(-14.9, -7.6, -0.3), SVector{3,Float64}(-22.2, -29.5, -36.8), false)]
yref1 = yref[1]
yref2 = yref[2]

@test isapprox(y1.cv, yref1.cv, atol=mctol)
@test isapprox(y1.cc, yref1.cc, atol=mctol)
@test isapprox(y1.Intv.lo, yref1.Intv.lo, atol=mctol)
@test isapprox(y1.Intv.hi, yref1.Intv.hi, atol=mctol)
@test isapprox(y1.cv_grad[1], yref1.cv_grad[1], atol=mctol)
@test isapprox(y1.cv_grad[2], yref1.cv_grad[2], atol=mctol)
@test isapprox(y1.cv_grad[3], yref1.cv_grad[3], atol=mctol)
@test isapprox(y1.cnst, yref1.cnst, atol=mctol)

@test isapprox(y2.cv, yref2.cv, atol=mctol)
@test isapprox(y2.cc, yref2.cc, atol=mctol)
@test isapprox(y2.Intv.lo, yref2.Intv.lo, atol=mctol)
@test isapprox(y2.Intv.hi, yref2.Intv.hi, atol=mctol)
@test isapprox(y2.cv_grad[1], yref2.cv_grad[1], atol=mctol)
@test isapprox(y2.cv_grad[2], yref2.cv_grad[2], atol=mctol)
@test isapprox(y2.cv_grad[3], yref2.cv_grad[3], atol=mctol)
@test isapprox(y2.cnst, yref2.cnst, atol=mctol)

end
