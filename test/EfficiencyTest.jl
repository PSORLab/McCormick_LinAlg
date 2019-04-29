module EfficiencyTest

    using Compat
    using Compat.Test
    using BenchmarkTools
    import BenchmarkTools.minimum
    using StaticArrays, SparseArrays
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

    include("../src/BLAS_MC/Lin_Functions/form.jl")
    m1 = MC{4}(4.0, 5.0, IntervalType(4,7), SVector{4,Float64}(3.0, 2.0, 4.4, 1.0), SVector{4,Float64}(3.0, 2.0, 4.0, 1.0), false)
    m2 = MC{4}(33.7,50.0,IntervalType(32,50),SVector{4,Float64}(64.0, 23.5, 8.0, 96.0),SVector{4,Float64}(54.0, 7.8, 3.6, 18.0),false)
    m3 = MC{4}(4.0, 5.0, IntervalType(4, 5), SVector{4,Float64}(4.0, 5.0, 13.2, 6.0), SVector{4,Float64}(3.0, 2.0, 5.2, 1.0),false)
    m4 = MC{4}(3.0, 4.0, IntervalType(2, 8), SVector{4,Float64}(4.0, 5.0, 7.0, 6.0), SVector{4,Float64}(3.0, 2.0, 2.5, 1.0), false)
    a = 6.3

    M = [m1,m2,m3,m4]
    Random.seed!(0);
    n = 75 #Vector Size for general
    ind = rand(1:4, n*2)
    X = map(x -> M[x], ind[1:n])
    Y = map(x -> M[x], ind[n+1:2n])

    ns = 1000 #Size for sparse systems

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
    println("GBMV efficiency") #Recently increased dimensions on this test
    ns =ms= 1000 #Size for sparse systems
    inds = rand(1:4, ns*2)
    AB = rand(M, ns,ns);
    xs = rand(M, ns);
    ys = rand(M, ns);
    kl = 100 #Efficiency is very dependent on the bandwidth
    ku = 100
    MCzero = zero(MC{4})
    for i in 1:ms #MAKE A a SPARSE BANDED HERE
        for j in 1:ns
            if j<i-kl || j>i+ku
                AB[i,j] = MCzero
            end
        end
    end

    MCrows = []
    MCcols = []
    ABs_ = []
    for i in 1:ns
           for j in 1:ns
           if AB[i,j] != MCzero
               push!(ABs_, AB[i,j])
               push!(MCrows, i)
               push!(MCcols, j)
           end
        end
    end
    ABs = sparse(MCcols, MCrows, ABs_)
    ABs = SparseMatrixCSC{MC{4}, Int64}(ABs)
    ABg = band(AB,ms,ns,ku,kl)

    bgbmvs = @benchmarkable GBMVs($TRANS, $alpha, $ABs, $xs, $beta, $ys)
    #bdsgbmvs = @benchmarkable deadsimplegbmvs($TRANS, $alpha, $ABs, $x, $beta, $y)
    bgbmv = @benchmarkable GBMV($TRANS, $ms, $ns, $kl, $ku, $alpha, $ABg, $xs, $beta, $ys)
    bdsgbmv = @benchmarkable deadsimplegemv($TRANS, $ms, $ns, $alpha, $AB, $xs, $beta, $ys) #gemv is simple matrix mult,
    for b in [bgbmv, bdsgbmv, bgbmvs]#, bdsgbmvs]                                                           #So it is benchmark for all special matrix forms
        tune!(b)
    end
    new, old, sp#=, oldsp =#= minimum(run(bgbmv)), minimum(run(bdsgbmv)), minimum(run(bgbmvs))#, minimum(run(bdsgbmvs))
    println(judge(new, old))
    println("SparseArray")
    println(judge(sp, old))
    #println(judge(sp, oldsp))

#SYMV
println("SYMV efficiency")
UPLO = "U"
AS = copy(A)
for i in 1:m #Make Symmetric (copy Upper to Lower)
    for j in 1:n
        if j < i
            AS[i,j] = AS[j,i]
        end
    end
end


bsymv = @benchmarkable SYMV($UPLO, $n, $alpha, $AS, $x, $beta, $y)
bdssymv = @benchmarkable deadsimplegemv($TRANS, $m, $n, $alpha, $AS, $x, $beta, $y) #gemv is simple matrix mult,
for b in [bsymv, bdssymv]                                                           #So it is benchmark for all special matrix forms
    tune!(b)
end
new, old = minimum(run(bsymv)), minimum(run(bdssymv))
println(judge(new, old))

#SBMV
println("SBMV efficiency")
DIAG = "N"
k = 100 #1000x1000 matrix and k=100 is 81% sparse
for i in 1:ms #Make AB Symmetric (copy Upper to Lower). Bandwidth is same
    for j in 1:ns
        if j < i
            AB[i,j] = AB[j,i]
        end
    end
end

MCrows = []
MCcols = []
ABs_ = []
for i in 1:ns
       for j in 1:ns
       if AB[i,j] != MCzero
           push!(ABs_, AB[i,j])
           push!(MCrows, i)
           push!(MCcols, j)
       end
    end
end
ABs = sparse(MCcols, MCrows, ABs_)
ABs = SparseMatrixCSC{MC{4}, Int64}(ABs)

ABu = sbandu(AB,ms,ns,k)

bgbmvs = @benchmarkable GBMVs($TRANS, $alpha, $ABs, $xs, $beta, $ys)
bsbmv = @benchmarkable SBMV($UPLO, $ns, $k, $alpha, $ABu, $xs, $beta, $ys)
bdssbmv = @benchmarkable deadsimplegemv($TRANS, $ms, $ns, $alpha, $AB, $xs, $beta, $ys) #gemv is simple matrix mult,
for b in [bsbmv, bdssbmv]                                                           #So it is benchmark for all special matrix forms
    tune!(b)
end
new, old , newsp= minimum(run(bsbmv)), minimum(run(bdssbmv)), minimum(run(bgbmvs))
println(judge(new, old))
println("Using Sparse Array")
println(judge(newsp, old))

#TRMV
println("TRMV efficiency")
DIAG = "N" #Whether matrix is unit ("U") triag or not
AT = copy(A)
for i in 1:m #Make Upper Triangular
    for j in 1:n
        if j < i
            AT[i,j] = MCzero
        end
    end
end
btrmv = @benchmarkable TRMV($UPLO, $TRANS, $DIAG, $n, $AT, $x)
bdstrmv = @benchmarkable deadsimplegemv($TRANS, $m, $n, $alpha, $AT, $x, $beta, $y) #gemv is simple matrix mult,
for b in [btrmv, bdstrmv]                                                           #So it is benchmark for all special matrix forms
    tune!(b)
end
new, old = minimum(run(btrmv)), minimum(run(bdstrmv))
println(judge(new, old))

end
