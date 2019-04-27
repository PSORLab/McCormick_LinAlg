module McCormick_LinAlg

    using EAGO
    export XSCAL, DOT, SAXPY, XNRM2, GEMV, GBMV, SYMV, SBMV, TRMV
    export testcorrectness
    export testefficiency

#Just include files individually for now, under this one module no REEXPORT =#

    using StaticArrays, SparseArrays
    include("BLAS_MC/Lin_Functions/DOT.jl")
    include("BLAS_MC/Lin_Functions/XSCAL.jl")
    include("BLAS_MC/Lin_Functions/XNRM2.jl")
    include("BLAS_MC/Lin_Functions/SAXPY.jl")
    include("BLAS_MC/Lin_Functions/GEMV.jl")
    include("BLAS_MC/Lin_Functions/GBMV.jl")
    include("BLAS_MC/Lin_Functions/SYMV.jl")
    include("BLAS_MC/Lin_Functions/SBMV.jl")
    include("BLAS_MC/Lin_Functions/TRMV.jl")
    source = @__DIR__
    function testcorrectness()              #This is calling from current working directory, not source of file
        include(source * "\\..\\test\\BLAS_MCtest.jl") #There was a trick to find package location in MathOptInt for this kind of call
    end
    function testefficiency()
        include(source * "\\..\\test\\EfficiencyTest.jl")
    end
end
