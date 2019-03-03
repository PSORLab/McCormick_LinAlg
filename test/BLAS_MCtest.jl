module Check_BLAS_MC

    using Compat
    using Compat.Test
    using StaticArrays
    using McCormick_LinAlg
    #using...

#    include("dotproduct.jl")
    include("xscal.jl")
end
