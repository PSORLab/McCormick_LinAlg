#@test lines for XSCAL etc
#Meant to be the easiest test file to make
@testset "Test XNRM2" begin

mctol = 1E-4
#= m = MC{2}(2.0, 3.0, IntervalType(1.0,4.0),
             seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false) =#
#To do make test vectors less similar
m1= MC{3}(4.0, 5.0, IntervalType(4.,5.), SVector{3,Float64}(4.0, 5.0, 13.0),SVector{3,Float64}(3.0, 12.5, 1.0), false)
m2= MC{3}(4.0, 5.0, IntervalType(6.,7.), SVector{3,Float64}(7.0, 2.0, 6.0), SVector{3,Float64}(6.0, 2.0, 1.0), false)
m3= MC{3}(4.0, 5.0, IntervalType(4., 5.),SVector{3,Float64}(4.0, 9.0, 6.0), SVector{3,Float64}(10.0, 7.0, 1.0),false)
m4= MC{3}(3.0, 4.0, IntervalType(3., 4.),SVector{3,Float64}(3.0, 5.0, 8.0), SVector{3,Float64}(3.0, 2.0, 2.0), false)
X = SVector{4,MC}(m1,m2,m3,m4)

yref = MC{3}(8.774964, 8.774964, IntervalType(8.77496, 10.7239), SVector{3,Float64}(4.20539, 7.28251, 10.2571), SVector{3,Float64}(4.67238, 8.0912, 11.3961), false)
#Lots of mixed results sqrting MC using sqrt and pow with .5
y = XNRM2(X)

@test isapprox(y.cv, yref.cv, atol=mctol)
@test isapprox(y.cc, yref.cc, atol=mctol)
@test isapprox(y.Intv.lo, yref.Intv.lo, atol=mctol)
@test isapprox(y.Intv.hi, yref.Intv.hi, atol=mctol)
@test isapprox(y.cv_grad[1], yref.cv_grad[1], atol=mctol)
@test isapprox(y.cv_grad[2], yref.cv_grad[2], atol=mctol)
@test isapprox(y.cv_grad[3], yref.cv_grad[3], atol=mctol)
@test isapprox(y.cnst, yref.cnst, atol=mctol)

end
