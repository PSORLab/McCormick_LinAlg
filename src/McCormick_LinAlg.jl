module McCormick_LinAlg

#=    using LinearAlgebra
    using SparseArrays

    using Reexport, StaticArrays

    import IntervalArithmetic: +, -, *, /, convert, in, isempty, one, zero,
                               real, eps, max, min, abs, exp,
                               expm1, log, log2, log10, log1p, sqrt,
                               sin, cos, tan, min, max, sec, csc, cot, step,
                               sign, dist, mid, pow, Interval, sinh, cosh, ∩,
                               IntervalBox, bisect, isdisjoint
=#
    #Load Package EAGO xEwhJ version commented out for just McCormick.jl
    using EAGO
    export XSCAL, DOT, SAXPY, XNRM2, GEMV
    export testcorrectness
    export testefficiency
#=    include("BLAS_MC/BLAS_MC.jl")
    using .BLAS_MC

Just include files individually for now, under this one module no REEXPORT =#
    #include("BLAS_MC/Lin_Functions/DOT.jl")
    using StaticArrays #should be redundant
    include("BLAS_MC/Lin_Functions/DOT.jl")
    include("BLAS_MC/Lin_Functions/XSCAL.jl")
    include("BLAS_MC/Lin_Functions/XNRM2.jl")
    include("BLAS_MC/Lin_Functions/SAXPY.jl")
    include("BLAS_MC/Lin_Functions/GEMV.jl")

    function testcorrectness()
        include("../test/BLAS_MCtest.jl")
    end
    function testefficiency()
        include("../test/EfficiencyTest")
    end

end
