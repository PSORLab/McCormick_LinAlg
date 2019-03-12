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
#=
With MC{N} where N=3

DOT efficiency
TrialJudgement(+107.76% => regression)
TrialJudgement(+899.00% => regression)
SAXPY efficiency
TrialJudgement(-40.74% => improvement)
TrialJudgement(+184.45% => regression)
XSCAL efficiency
TrialJudgement(+8.57% => regression)
TrialJudgement(+1372.66% => regression)

With N = 50
DOT efficiency
TrialJudgement(+41.81% => regression)
TrialJudgement(+796.48% => regression)
SAXPY efficiency
TrialJudgement(+197.74% => regression)
TrialJudgement(+1215.94% => regression)
XSCAL efficiency
TrialJudgement(+2328.73% => regression)
TrialJudgement(+3137.82% => regression)
=#
    #using...
#Efficient Functions

#Simple Implementations of each function
    include("SimpleFunctions/simpledot.jl")
    include("SimpleFunctions/simplexscal.jl")
    include("SimpleFunctions/simplesaxpy.jl")
    include("SimpleFunctions/simplexnrm2.jl")

    G = BenchmarkGroup()
    G["opt"] = BenchmarkGroup(["optimized", "BLAS"])
    G["bench"] = BenchmarkGroup(["ineffifient", "slow", "simple"])

    m1 = MC{3}(4.0, 5.0, IntervalType(4,7), SVector{3,Float64}(3.0, 2.0, 1.0), SVector{3,Float64}(3.0, 2.0, 1.0), false)
    m2 = MC{3}(3.2,50.0,IntervalType(32,50),SVector{3,Float64}(64.0,8.0, 96.0),SVector{3,Float64}(54.0,3.6, 18.0),false)
    m3 = MC{3}(4.0, 5.0, IntervalType(4, 5), SVector{3,Float64}(4.0, 5.0, 6.0), SVector{3,Float64}(3.0, 2.0, 1.0),false)
    m4 = MC{3}(3.0, 4.0, IntervalType(3, 4), SVector{3,Float64}(4.0, 5.0, 6.0), SVector{3,Float64}(3.0, 2.0, 1.0), false)
    a = 6.3

    M = (m1,m2,m3,m4)
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

    bxscal = @benchmarkable XSCAL($X, $a)
    bsxscal = @benchmarkable simplexscal($X, $a)
    bdsxscal = @benchmarkable deadsimplexscal($X, $a)
    for b in [bxscal, bsxscal, bdsxscal]
        tune!(b)
    end
    new, old1, old2 = minimum(run(bxscal)), minimum(run(bsxscal)), minimum(run(bdsxscal))
    println(judge(new, old1))
    println(judge(new, old2))

end
