#NOT COMPLETE
#Symmetric Banded Matrix Vector Multiply
@testset "Test SBMV" begin

mctol = 2E-3
#Since errors are rounding, may just need to enter more SigFigs in test data
 m1 = MC{3}(5.0, 13.0, IntervalType(4,15), SVector{3,Float64}([4.0, 3.0, 18.0]), SVector{3,Float64}([11.0, 12.0, 8.0]), false)
 m2 = MC{3}(2.0, 3.0, IntervalType(1,7), SVector{3,Float64}([2.0, 16.0, 17.0]), SVector{3,Float64}([12.0, 13.0, 11.0]), false)
 m3 = MC{3}(8.0, 16.0, IntervalType(6,20), SVector{3,Float64}([16.0, 16.0, 4.0]), SVector{3,Float64}([15.0, 7.0, 10.0]), false)
 m4 = MC{3}(9.0, 15.0, IntervalType(3,18), SVector{3,Float64}([10.0, 8.0, 4.0]), SVector{3,Float64}([16.0, 9.0, 13.0]), false)
#Test R^8x10 * R^10 + R^8
#2 lower- and 2 upper-diagonals
M = [m1,m2,m3,m4]
Random.seed!(0)

m = 10
n = 10
k = 2
AF = rand(M, m,n)
MCzero = MC{3}(0.0,0.0)
for i in 1:m #MAKE A a symmetric SPARSE BANDED HERE
    for j in 1:n
        if j<i-k || j>i+k
            AF[i,j] = MCzero
        end
    end
end
for i in 1:m #Make Symmetric (copy Upper to Lower)
    for j in 1:n
        if j < i
            AF[i,j] = AF[j,i]
        end
    end
end
include("../src/BLAS_MC/Lin_Functions/form.jl")
A = sbandu(AF,m,n,k) #AF is reformatted to a banded matrix storage form fit for this function

x = rand(M, n)
y_ = rand(M, m)
alpha, beta = 2.0, 6.1
UPLO = "U"
#y_ref = alpha*AF*x + beta*y_;
y = SBMV(UPLO, n, k, alpha, A, x, beta, y_)


testset = [1,3,6,10]
y1, y2, y3, y4 = map(i -> y[i], testset)

yref1=MC{3}(256.5, 1169.3, IntervalType(96.3999, 1631.5), SVector{3,Float64}(400.4, 516.3, 533.8), SVector{3,Float64}(1855.1, 1309.2, 1462.8), false)
yref2=MC{3}(298.5, 1357.3, IntervalType(110.399, 2093.5), SVector{3,Float64}(456.4, 762.3, 815.8), SVector{3,Float64}(2701.1, 2209.2, 2230.8), false)
yref3=MC{3}(160.2, 806.3, IntervalType(48.0999, 1376.71), SVector{3,Float64}(228.2, 679.6, 671.7), SVector{3,Float64}(2021.2, 1821.3, 1771.1), false)
yref4=MC{3}(208.9, 829.5, IntervalType(104.299, 1287.81), SVector{3,Float64}(509.0, 720.8, 400.4), SVector{3,Float64}(1579.6, 1160.9, 1235.3), false)



@test isapprox(y1.Intv.lo, yref1.Intv.lo, atol = mctol)
@test isapprox(y1.Intv.hi, yref1.Intv.hi, atol = mctol)

@test isapprox(y2.Intv.lo, yref2.Intv.lo, atol = mctol)
@test isapprox(y2.Intv.hi, yref2.Intv.hi, atol = mctol)

@test isapprox(y3.Intv.lo, yref3.Intv.lo, atol = mctol)
@test isapprox(y3.Intv.hi, yref3.Intv.hi, atol = mctol)

@test isapprox(y4.Intv.lo, yref4.Intv.lo, atol = mctol)
@test isapprox(y4.Intv.hi, yref4.Intv.hi, atol = mctol)

UPLO = "L"
A = sbandl(AF,m,n,k) #AF is reformatted to a banded matrix storage form fit for this function

y = SBMV(UPLO, n, k, alpha, A, x, beta, y_)#Using the lower triangular of A

        testset = [1,3,6,10]
        y1, y2, y3, y4 = map(i -> y[i], testset)

#y_ref are the same

@test isapprox(y1.Intv.lo, yref1.Intv.lo, atol = mctol)
@test isapprox(y1.Intv.hi, yref1.Intv.hi, atol = mctol)

@test isapprox(y2.Intv.lo, yref2.Intv.lo, atol = mctol)
@test isapprox(y2.Intv.hi, yref2.Intv.hi, atol = mctol)

@test isapprox(y3.Intv.lo, yref3.Intv.lo, atol = mctol)
@test isapprox(y3.Intv.hi, yref3.Intv.hi, atol = mctol)

@test isapprox(y4.Intv.lo, yref4.Intv.lo, atol = mctol)
@test isapprox(y4.Intv.hi, yref4.Intv.hi, atol = mctol)

end
