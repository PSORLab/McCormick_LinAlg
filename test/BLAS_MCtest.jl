module Check_BLAS_MC

    using Compat
    using Compat.Test
    using StaticArrays
    using EAGO #Should be dependent in LinAlg So idk why MC not defined
    using McCormick_LinAlg
    using Random

#BLAS Level 1
    include("dot.jl")
    include("xscal.jl")
    include("saxpy.jl")
    include("xnrm2.jl")
#BLAS Level 2
    include("gemv.jl")
    include("gbmv.jl")
    include("symv.jl")
    include("sbmv.jl")
    include("trmv.jl")
    #include("tbmv.jl") #Not made yet
#BLAS Level 3
end
