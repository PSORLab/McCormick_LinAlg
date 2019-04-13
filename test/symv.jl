#NOT COMPLETE - Need formatting for matrix, test data

@testset "Test SYMV" begin

mctol = 2E-3

m1 = MC{3}(5.0, 13.0, IntervalType(4,15), SVector{3,Float64}([4.0, 3.0, 18.0]), SVector{3,Float64}([11.0, 12.0, 8.0]), false)
m2 = MC{3}(2.0, 3.0, IntervalType(1,7), SVector{3,Float64}([2.0, 16.0, 17.0]), SVector{3,Float64}([12.0, 13.0, 11.0]), false)
m3 = MC{3}(8.0, 16.0, IntervalType(6,20), SVector{3,Float64}([16.0, 16.0, 4.0]), SVector{3,Float64}([15.0, 7.0, 10.0]), false)
m4 = MC{3}(9.0, 15.0, IntervalType(3,18), SVector{3,Float64}([10.0, 8.0, 4.0]), SVector{3,Float64}([16.0, 9.0, 13.0]), false)

M = [m1,m2,m3,m4]
Random.seed!(0)

m = 10
n = 10
AF = rand(M, m,n)
MCzero = MC{3}(0.0,0.0)
for i in 1:m #Make Symmetric (copy Upper to Lower)
    for j in 1:n
        if j < i
            AF[i,j] = AF[j,i]
        end
    end
end
x = rand(M, n)
y_ = rand(M, m)
alpha, beta = 2.0, 6.1
UPLO = "U"

y = SYMV(UPLO, n, alpha, AF, x, beta, y_)
testset = [1,3,6,10]
y1, y2, y3, y4 = map(i -> y[i], testset)
y_ref = alpha*AF*x + beta*y_
yref1, yref2, yref3, yref4 = map(i -> y_ref[i], testset)

@test isapprox(y1.Intv.lo, yref1.Intv.lo, atol = mctol)
@test isapprox(y2.Intv.hi, yref2.Intv.hi, atol = mctol)
@test isapprox(y3.Intv.lo, yref3.Intv.lo, atol = mctol)
@test isapprox(y4.Intv.hi, yref4.Intv.hi, atol = mctol)
###################################################################
UPLO = "L" #Using lower triangular of A for calulcations

y = SYMV(UPLO, n, alpha, AF, x, beta, y_)
testset = [1,3,6,10]
y1, y2, y3, y4 = map(i -> y[i], testset)
yref1, yref2, yref3, yref4 = map(i -> y_ref[i], testset)

@test isapprox(y1.Intv.lo, yref1.Intv.lo, atol = mctol)
@test isapprox(y2.Intv.hi, yref2.Intv.hi, atol = mctol)
@test isapprox(y3.Intv.lo, yref3.Intv.lo, atol = mctol)
@test isapprox(y4.Intv.hi, yref4.Intv.hi, atol = mctol)
end
