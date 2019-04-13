module EfficiencyTest

    using Compat
    using Compat.Test
    using BenchmarkTools
    import BenchmarkTools.minimum
    using StaticArrays
    using EAGO #Should be dependent in LinAlg So idk why MC not defined
    using McCormick_LinAlg
    using Random
    #import McCormick_LinAlg.XSCAL #shouldnt be necessary with proper export

    #using...
#Efficient Functions
# Can still implement @code_warntype, code_llvm, @inbounds, @fastmath, @simd,
#Consider whether recycling functions into eachother will speed up or not.
#Simple Implementations of each function
    include("SimpleFunctions/simpledot.jl")
    include("SimpleFunctions/simplexscal.jl")
    include("SimpleFunctions/simplesaxpy.jl")
    include("SimpleFunctions/simplexnrm2.jl")
    include("SimpleFunctions/simplegemv.jl")

    m1 = MC{3}(4.0, 5.0, IntervalType(4,7), SVector{3,Float64}(3.0, 2.0, 1.0), SVector{3,Float64}(3.0, 2.0, 1.0), false)
    m2 = MC{3}(33.7,50.0,IntervalType(32,50),SVector{3,Float64}(64.0,8.0, 96.0),SVector{3,Float64}(54.0,3.6, 18.0),false)
    m3 = MC{3}(4.0, 5.0, IntervalType(4, 5), SVector{3,Float64}(4.0, 5.0, 6.0), SVector{3,Float64}(3.0, 2.0, 1.0),false)
    m4 = MC{3}(3.0, 4.0, IntervalType(2, 8), SVector{3,Float64}(4.0, 5.0, 6.0), SVector{3,Float64}(3.0, 2.0, 1.0), false)
    a = 6.3

    M = [m1,m2,m3,m4]
    Random.seed!(0)
    n = 50 #Vector Size for all testing
    ind = rand(1:4, n*2)
    X = map(x -> M[x], ind[1:n])
    Y = map(x -> M[x], ind[n+1:2n])


#DOT
    println("DOT efficiency")

    bdot = @benchmarkable DOT($X, $Y)
    bsdot = @benchmarkable simpledot($X, $Y)
    bdsdot = @benchmarkable deadsimpledot($X, $Y)
    for b in [bdot, bsdot, bdsdot]
        tune!(b)
    end
    new, old1, old2 = minimum(run(bdot)), minimum(run(bsdot)), minimum(run(bdsdot))
    println(judge(new, old1))
    println(judge(new, old2))

#SAXPY
    println("SAXPY efficiency")

    bsaxpy = @benchmarkable SAXPY($a, $X, $Y)
    bssaxpy = @benchmarkable simplesaxpy($a, $X, $Y)
    bdssaxpy = @benchmarkable deadsimplesaxpy($a, $X, $Y)
    for b in [bsaxpy, bssaxpy, bdssaxpy]
        tune!(b)
    end
    new, old1, old2 = minimum(run(bsaxpy)), minimum(run(bssaxpy)), minimum(run(bdssaxpy))
    println(judge(new, old1))
    println(judge(new, old2))

#XSCAL
    println("XSCAL efficiency")

    bxscal = @benchmarkable XSCAL($a, $X)
    bsxscal = @benchmarkable simplexscal($a, $X)
    bdsxscal = @benchmarkable deadsimplexscal($a, $X)
    for b in [bxscal, bsxscal, bdsxscal]
        tune!(b)
    end
    new, old1, old2 = minimum(run(bxscal)), minimum(run(bsxscal)), minimum(run(bdsxscal))
    println(judge(new, old1))
    println(judge(new, old2))

    #GEMV
    println("GEMV efficiency")

    M = [m1,m2,m3,m4];
    m = n;
    A = rand(M, m,n);
    x = rand(M, n);
    y = rand(M, m);
    alpha, beta = 2.0, 6.1;
    TRANS = "N";

    bgemv = @benchmarkable GEMV($TRANS, $m, $n, $alpha, $A, $x, $beta, $y)
    bdsgemv = @benchmarkable deadsimplegemv($TRANS, $m, $n, $alpha, $A, $x, $beta, $y)
    for b in [bgemv, bdsgemv]
        tune!(b)
    end
    new, old2 = minimum(run(bgemv)), minimum(run(bdsgemv))
    println(judge(new, old2))

#GBVM
    AB = copy(A)
    kl = 5 #Efficiency is very dependednt on the bandwidth
    ku = 4
    AB = rand(M, m,n)
    MCzero = MC{3}(0.0,0.0)
    for i in 1:m #MAKE A a SPARSE BANDED HERE
        for j in 1:n
            if j<i-kl || j>i+ku
                AB[i,j] = MCzero
            end
        end
    end
    include("../../src/BLAS_MC/Lin_Functions/form.jl")
    ABg = band(AB,m,n,ku,kl)
    bgbmv = @benchmarkable GBMV($TRANS, $m, $n, $kl, $ku, $alpha, $ABg, $x, $beta, $y)
    bdsgbmv = @benchmarkable deadsimplegemv($TRANS, $m, $n, $alpha, $AB, $x, $beta, $y) #gemv is simple matrix mult,
    for b in [bgbmv, bdsgbmv]                                                           #So it is benchmark for all special matrix forms
        tune!(b)
    end
    new, old = minimum(run(bgmbv)), minimum(run(bdsgbmv))
    println(judge(new, old))


end
