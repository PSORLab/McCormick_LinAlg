module Check_BLAS_MC

    using Compat
    using Compat.Test
    using StaticArrays
    using EAGO #Should be dependent in LinAlg So idk why MC not defined
    using McCormick_LinAlg
    using Random

#BLAS Level 1
    include("BLAS_MC/dot.jl")
    include("BLAS_MC/xscal.jl")
    include("BLAS_MC/saxpy.jl")
    include("BLAS_MC/xnrm2.jl")
#BLAS Level 2
    include("BLAS_MC/gemv.jl")
    include("BLAS_MC/gbmv.jl")
    include("BLAS_MC/symv.jl")
    include("BLAS_MC/sbmv.jl")
    include("BLAS_MC/trmv.jl")
    #include("tbmv.jl") #Not made yet
#BLAS Level 3
end
