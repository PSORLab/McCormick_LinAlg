#Need 5+ size vector to test
#Meant to be the easiest test file to make
@testset "Test DOT" begin
mctol = 1E-4
#= m = MC{2}(2.0, 3.0, IntervalType(1.0,4.0),
             seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false) =#
m1 = MC{3}(4.0, 5.0, IntervalType(4,4), SVector{3,Float64}(3.0, 2.0, 1.0), SVector{3,Float64}(3.0, 2.0, 1.0), false)
m2 = MC{3}(3.2,50.0,IntervalType(32,50),SVector{3,Float64}(64.0,8.0, 96.0),SVector{3,Float64}(54.0,3.6, 18.0),false)
mv1 = SVector{2,MC}(m1,m2)
m3 = MC{3}(4.0, 5.0, IntervalType(4, 5), SVector{3,Float64}(4.0, 5.0, 6.0), SVector{3,Float64}(3.0, 2.0, 1.0),false)
m4 = MC{3}(3.0, 4.0, IntervalType(3, 4), SVector{3,Float64}(4.0, 5.0, 6.0), SVector{3,Float64}(3.0, 2.0, 1.0), false)
mv2 = SVector{2,MC}(m3,m4)

y = DOT(mv1, mv2)
yref = MC{3}(112.0, 220.0, IntervalType(112, 220), SVector{3,Float64}(28.0, 28.0, 28.0), SVector{3,Float64}(312.0, 78.4, 104.0), false)

@test y.cv == yref.cv
@test y.cc == yref.cc
@test y.Intv.lo == yref.Intv.lo
@test y.Intv.hi == yref.Intv.hi
@test y.cv_grad[1] == yref.cv_grad[1]
@test y.cv_grad[2] == yref.cv_grad[2]
@test y.cv_grad[3] == yref.cv_grad[3]
@test y.cnst == yref.cnst

m1 = MC{3}(4.0, 5.0, IntervalType(4,4), SVector{3,Float64}(3.0, 2.0, 1.0), SVector{3,Float64}(3.0, 2.0, 1.0), false)
m2 = MC{3}(3.2,50.0,IntervalType(32,50),SVector{3,Float64}(64.0,8.0, 96.0),SVector{3,Float64}(54.0,3.6, 18.0),false)
mv1 = SVector{12, MC}(m1,m1,m1,m1, m1,m1,m1,m1, m1,m1,m1,m1)
mv2 = SVector{12, MC}(m2,m2,m2,m2, m2,m2,m2,m2, m2,m2,m2,m2)

y = DOT(mv1,mv2)
yref = MC{3}(1536.0, 2400.0, IntervalType(1536, 2400), SVector{3,Float64}(0.0, 0.0, 0.0), SVector{3,Float64}(0.0, 0.0, 0.0), false)

@test y.cv == yref.cv
@test y.cc == yref.cc
@test y.Intv.lo == yref.Intv.lo
@test y.Intv.hi == yref.Intv.hi
@test y.cv_grad[1] == yref.cv_grad[1]
@test y.cv_grad[2] == yref.cv_grad[2]
@test y.cv_grad[3] == yref.cv_grad[3]
@test y.cnst == yref.cnst

end
