module Check_BLAS_MC

    using Compat
    using Compat.Test

    using McCormick_LinAlg, StaticArrays, IntervalArithmetic
    #using...

    include("scalar_routines.jl")
    include("dotproduct.jl")
end
