module EfficiencyTest

    using Compat
    using Compat.Test
    using BenchmarkTools
    import BenchmarkTools.mean
    using StaticArrays
    using EAGO #Should be dependent in LinAlg So idk why MC not defined
    using McCormick_LinAlg
    #import McCormick_LinAlg.XSCAL #shouldnt be necessary with proper export

    #using...
#Efficient Functions

#Simple Implementations of each function
    include("SimpleFunction/simpledot.jl")
    include("SimpleFunction/simplexscal.jl")
    include("SimpleFunction/simplesaxpy.jl")
    include("SimpleFunction/simplexnrm2.jl")

    G = BenchmarkGroup()
    G["opt"] = BenchmarkGroup(["optimized", "BLAS"])
    G["bench"] = BenchmarkGroup(["ineffifient", "slow", "simple"])
#DOT
    println("DOT efficiency")
    m1 = MC{3}(4.0, 5.0, IntervalType(4,4), SVector{3,Float64}(3.0, 2.0, 1.0), SVector{3,Float64}(3.0, 2.0, 1.0), false)
    m2 = MC{3}(3.2,50.0,IntervalType(32,50),SVector{3,Float64}(64.0,8.0, 96.0),SVector{3,Float64}(54.0,3.6, 18.0),false)
    mv1 = SVector{2,MC}(m1,m2)
    m3 = MC{3}(4.0, 5.0, IntervalType(4, 5), SVector{3,Float64}(4.0, 5.0, 6.0), SVector{3,Float64}(3.0, 2.0, 1.0),false)
    m4 = MC{3}(3.0, 4.0, IntervalType(3, 4), SVector{3,Float64}(4.0, 5.0, 6.0), SVector{3,Float64}(3.0, 2.0, 1.0), false)
    mv2 = SVector{2,MC}(m3,m4)

    bdot = @benchmarkable DOT($mv1, $mv2)
    bsdot = @benchmarkable simpledot($mv1, $mv2)
    bdsdot = @benchmarkable deadsimpledot($mv1, $mv2)
    for b in [bdot, bsdot, bdsdot]
        tune!(b)
    end
    [new, old1, old2] = [run(bdot), run(bsdot), run(bdsdot)]
    judge(new, old1)
    judge(new, old2)

#SAXPY
    println("SAXPY efficiency")
    m1= MC{3}(4.0, 5.0, IntervalType(4.,5.), SVector{3,Float64}(4.0, 5.0, 6.0), SVector{3,Float64}(3.0, 2.0, 1.0), false)
    m2= MC{3}(4.0, 5.0, IntervalType(4.,5.), SVector{3,Float64}(4.0, 5.0, 6.0), SVector{3,Float64}(3.0, 2.0, 1.0), false)
    m3= MC{3}(4.0, 5.0, IntervalType(4., 5.), SVector{3,Float64}(4.0, 5.0, 6.0), SVector{3,Float64}(3.0, 2.0, 1.0),false)
    m4= MC{3}(3.0, 4.0, IntervalType(3., 4.), SVector{3,Float64}(4.0, 5.0, 6.0), SVector{3,Float64}(3.0, 2.0, 1.0), false)
    X = SVector{2,MC}(m1,m2)
    Y = SVector{2,MC}(m3, m4)
    a::Float64 = 6.3

    bsaxpy = @benchmarkable SAXPY(a, X, Y)
    bssaxpy = @benchmarkable simplesaxpy(a, X, Y)
    bdssaxpy = @benchmarkable deadsimplesaxpy(a, X, Y)
    for b in [bsaxpy, bssaxpy, bdssaxpy]
        tune!(b)
    end
    [new, old1, old2] = [run(bsaxpy), run(bssaxpy), run(bdssaxpy)]
    judge(new, old1)
    judge(new, old2)

#XSCAL
    println("XSCAL efficiency")
    X = MC{3}(4.0, 5.0, IntervalType(4.,5.), SVector{3,Float64}(4.0, 5.0, 6.0), SVector{3,Float64}(3.0, 2.0, 1.0), false)
    a::Float64 = 6.3

    bxscal = @benchmarkable XSCAL(a, X)
    bsxscal = @benchmarkable simplexscal(a, X)
    bdsxscal = @benchmarkable deadsimplexscal(a, X)
    for b in [bxscal, bsxscal, bdsxscal]
        tune!(b)
    end
    [new, old1, old2] = [run(bxscal), run(bsxscal), run(bdsxscal)]
    judge(new, old1)
    judge(new, old2)

end
