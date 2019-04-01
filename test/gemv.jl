#@test lines for GEMV etc
#Rounding errors on Transpose case. All failures within 10%
#       "T" cv and cv_grad fields are wrong. Only 3rd out of 5 (middle) is completely right
@testset "Test GEMV" begin

mctol = 2E-3
#= m = MC{2}(2.0, 3.0, IntervalType(1.0,4.0),
             seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false) =#
m1 = MC{3}(5.0, 13.0, IntervalType(4,5), SVector{3,Float64}([4.0, 5.0, 6.0]), SVector{3,Float64}([3.0, 2.0, 1.0]), false)
m2 = MC{3}(1.0, 3.0, IntervalType(6,7), SVector{3,Float64}([12.0, 4.0, 8.0]), SVector{3,Float64}([2.0, 3.0, 4.0]), false)
m3 = MC{3}(4.0, 6.0, IntervalType(6,10), SVector{3,Float64}([1.0, 3.0, 8.0]), SVector{3,Float64}([2.0, 10.0, 4.0]), false)
m4 = MC{3}(6.0, 10.0, IntervalType(3,7), SVector{3,Float64}([12.0, 4.0, 3.0]), SVector{3,Float64}([1.0, 3.0, 10.0]), false)

#Test R^5x10 * R^10 + R^5

M = [m1,m2,m3,m4]
Random.seed!(0)

m = 10
n = 5
A = rand(M, m,n)
x = rand(M, n)
yparam = rand(M, m)
yc = copy(yparam)
alpha, beta = 2.0, 6.1
TRANS = "N"

y = GEMV(TRANS, m, n, alpha, A, x, beta, yparam) #Cheesey test, make unique solutions

yref = [MC{3}(252.1, 406.3, IntervalType(240.599, 450.701), SVector{3,Float64}([465.2, 206.4, 246.8]), SVector{3,Float64}([52.2, 96.3, 200.4]), false),
        MC{3}(306.1, 312.3, IntervalType(300.599, 592.701), SVector{3,Float64}([297.2, 374.4, 258.8]), SVector{3,Float64}([140.2, 518.3, 280.4]), false),
        MC{3}(282.4, 366.6, IntervalType(288.599, 539), SVector{3,Float64}([222.1, 90.3, 132.8]), SVector{3,Float64}([116.2, 319.0, 328.4]), false),
        MC{3}(348.6, 297.0, IntervalType(330.299, 662.701), SVector{3,Float64}([73.2, 24.4, 18.3]), SVector{3,Float64}([198.1, 558.3, 445.0]), false),
        MC{3}(258.1, 414.3, IntervalType(204.599, 478.701), SVector{3,Float64}([857.2, 388.4, 444.8]), SVector{3,Float64}([92.2, 174.3, 376.4]), false),
        MC{3}(324.6, 367.0, IntervalType(306.299, 580.701), SVector{3,Float64}([73.2, 24.4, 18.3]), SVector{3,Float64}([122.1, 416.3, 293.0]), false),
        MC{3}(282.6, 479.0, IntervalType(222.299, 480.701), SVector{3,Float64}([465.2, 206.4, 216.3]), SVector{3,Float64}([46.1, 96.3, 237.0]), false),
        MC{3}(318.6, 411.0, IntervalType(234.299, 560.701), SVector{3,Float64}([615.2, 272.4, 300.3]), SVector{3,Float64}([150.1, 452.3, 541.0]), false),
        MC{3}(278.1, 360.3, IntervalType(260.599, 528.701), SVector{3,Float64}([519.2, 262.4, 336.8]), SVector{3,Float64}([216.2, 440.3, 588.4]), false),
        MC{3}(268.1, 308.3, IntervalType(256.599, 480.701), SVector{3,Float64}([449.2, 296.4, 372.8]), SVector{3,Float64}([260.2, 456.3, 544.4]), false)]





testset = [1,3,6,10]

y1, y2, y3, y4 = map(i -> y[i], testset)
yref1, yref2, yref3, yref4 = map(i -> yref[i], testset)



@test isapprox(y1.cv, yref1.cv, atol = mctol)
@test isapprox(y1.cc, yref1.cc, atol = mctol)
@test isapprox(y1.Intv.lo, yref1.Intv.lo, atol = mctol)
@test isapprox(y1.Intv.hi, yref1.Intv.hi, atol = mctol)
@test isapprox(y1.cv_grad[1], yref1.cv_grad[1], atol = mctol)
@test isapprox(y1.cv_grad[2], yref1.cv_grad[2], atol = mctol)
@test isapprox(y1.cv_grad[3], yref1.cv_grad[3], atol = mctol)
@test isapprox(y1.cc_grad[1], yref1.cc_grad[1], atol = mctol)
@test isapprox(y1.cc_grad[2], yref1.cc_grad[2], atol = mctol)
@test isapprox(y1.cc_grad[3], yref1.cc_grad[3], atol = mctol)
@test isapprox(y1.cnst, yref1.cnst, atol = mctol)

@test isapprox(y2.cv, yref2.cv, atol = mctol)
@test isapprox(y2.cc, yref2.cc, atol = mctol)
@test isapprox(y2.Intv.lo, yref2.Intv.lo, atol = mctol)
@test isapprox(y2.Intv.hi, yref2.Intv.hi, atol = mctol)
@test isapprox(y2.cv_grad[1], yref2.cv_grad[1], atol = mctol)
@test isapprox(y2.cv_grad[2], yref2.cv_grad[2], atol = mctol)
@test isapprox(y2.cv_grad[3], yref2.cv_grad[3], atol = mctol)
@test isapprox(y2.cc_grad[1], yref2.cc_grad[1], atol = mctol)
@test isapprox(y2.cc_grad[2], yref2.cc_grad[2], atol = mctol)
@test isapprox(y2.cc_grad[3], yref2.cc_grad[3], atol = mctol)
@test isapprox(y2.cnst, yref2.cnst, atol = mctol)


@test isapprox(y3.cv, yref3.cv, atol = mctol)
@test isapprox(y3.cc, yref3.cc, atol = mctol)
@test isapprox(y3.Intv.lo, yref3.Intv.lo, atol = mctol)
@test isapprox(y3.Intv.hi, yref3.Intv.hi, atol = mctol)
@test isapprox(y3.cv_grad[1], yref3.cv_grad[1], atol = mctol) #These 3 all failed
@test isapprox(y3.cv_grad[2], yref3.cv_grad[2], atol = mctol)
@test isapprox(y3.cv_grad[3], yref3.cv_grad[3], atol = mctol)
@test isapprox(y3.cc_grad[1], yref3.cc_grad[1], atol = mctol)
@test isapprox(y3.cc_grad[2], yref3.cc_grad[2], atol = mctol)
@test isapprox(y3.cc_grad[3], yref3.cc_grad[3], atol = mctol)
@test isapprox(y3.cnst, yref3.cnst, atol = mctol)

@test isapprox(y4.cv, yref4.cv, atol = mctol)
@test isapprox(y4.cc, yref4.cc, atol = mctol)
@test isapprox(y4.Intv.lo, yref4.Intv.lo, atol = mctol)
@test isapprox(y4.Intv.hi, yref4.Intv.hi, atol = mctol)
@test isapprox(y4.cv_grad[1], yref4.cv_grad[1], atol = mctol)
@test isapprox(y4.cv_grad[2], yref4.cv_grad[2], atol = mctol)
@test isapprox(y4.cv_grad[3], yref4.cv_grad[3], atol = mctol)
@test isapprox(y4.cc_grad[1], yref4.cc_grad[1], atol = mctol)#These 3 failed
@test isapprox(y4.cc_grad[2], yref4.cc_grad[2], atol = mctol)
@test isapprox(y4.cc_grad[3], yref4.cc_grad[3], atol = mctol)
@test isapprox(y4.cnst, yref4.cnst, atol = mctol)

TRANS = "T"
yparam = yc #GEMV set!'s y so need to reset value (fixed)
x, yparam = yparam, x
y = GEMV(TRANS, m, n, alpha, A, x, beta, yparam)
#not communitive xy =/= yx for MC
#sometimes its bounds that work but arent the same
#both cv and cc should be within Interval bounds
#interval bounds may be a generally good thing to check
yref = [MC{3}(586.4, 714.6, IntervalType(474.599, 1055), SVector{3,Float64}([1100.1, 450.3, 498.8]), SVector{3,Float64}([260.2, 753.0, 808.4]), false),
        MC{3}(620.5, 777.3, IntervalType(468.399, 1168.5), SVector{3,Float64}([1494.4, 536.5, 540.6]), SVector{3,Float64}([358.3, 1000.2, 1070.1]), false),
        MC{3}(576.5, 819.3, IntervalType(480.399, 1122.5), SVector{3,Float64}([1122.4, 444.5, 540.6]), SVector{3,Float64}([374.3, 1032.2, 1294.1]), false),
        MC{3}(574.1, 734.3, IntervalType(438.599, 1094.71), SVector{3,Float64}([1559.2, 638.4, 696.8]), SVector{3,Float64}([300.2, 788.3, 984.4]), false),
        MC{3}(550.1, 844.3, IntervalType(438.599, 1092.71), SVector{3,Float64}([1507.2, 550.4, 636.8]), SVector{3,Float64}([304.2, 858.3, 1184.4]), false)]

        testset = [1,3,5]

        y1, y2, y3 = map(i -> y[i], testset)
        yref1, yref2, yref3 = map(i -> yref[i], testset)

@test isapprox(y1.cv, yref1.cv, atol = mctol)
@test isapprox(y1.cc, yref1.cc, atol = mctol)
@test isapprox(y1.Intv.lo, yref1.Intv.lo, atol = mctol)
@test isapprox(y1.Intv.hi, yref1.Intv.hi, atol = mctol)
@test isapprox(y1.cv_grad[1], yref1.cv_grad[1], atol = mctol)
@test isapprox(y1.cv_grad[2], yref1.cv_grad[2], atol = mctol)
@test isapprox(y1.cv_grad[3], yref1.cv_grad[3], atol = mctol)
@test isapprox(y1.cc_grad[1], yref1.cc_grad[1], atol = mctol)
@test isapprox(y1.cc_grad[2], yref1.cc_grad[2], atol = mctol)
@test isapprox(y1.cc_grad[3], yref1.cc_grad[3], atol = mctol)
@test isapprox(y1.cnst, yref1.cnst, atol = mctol)

@test isapprox(y2.cv, yref2.cv, atol = mctol)
@test isapprox(y2.cc, yref2.cc, atol = mctol)
@test isapprox(y2.Intv.lo, yref2.Intv.lo, atol = mctol)
@test isapprox(y2.Intv.hi, yref2.Intv.hi, atol = mctol)
@test isapprox(y2.cv_grad[1], yref2.cv_grad[1], atol = mctol)
@test isapprox(y2.cv_grad[2], yref2.cv_grad[2], atol = mctol)
@test isapprox(y2.cv_grad[3], yref2.cv_grad[3], atol = mctol)
@test isapprox(y2.cc_grad[1], yref2.cc_grad[1], atol = mctol)
@test isapprox(y2.cc_grad[2], yref2.cc_grad[2], atol = mctol)
@test isapprox(y2.cc_grad[3], yref2.cc_grad[3], atol = mctol)
@test isapprox(y2.cnst, yref2.cnst, atol = mctol)


@test isapprox(y3.cv, yref3.cv, atol = mctol)
@test isapprox(y3.cc, yref3.cc, atol = mctol)
@test isapprox(y3.Intv.lo, yref3.Intv.lo, atol = mctol)
@test isapprox(y3.Intv.hi, yref3.Intv.hi, atol = mctol)
@test isapprox(y3.cv_grad[1], yref3.cv_grad[1], atol = mctol) #These 3 all failed
@test isapprox(y3.cv_grad[2], yref3.cv_grad[2], atol = mctol)
@test isapprox(y3.cv_grad[3], yref3.cv_grad[3], atol = mctol)
@test isapprox(y3.cc_grad[1], yref3.cc_grad[1], atol = mctol)
@test isapprox(y3.cc_grad[2], yref3.cc_grad[2], atol = mctol)
@test isapprox(y3.cc_grad[3], yref3.cc_grad[3], atol = mctol)
@test isapprox(y3.cnst, yref3.cnst, atol = mctol)

end
