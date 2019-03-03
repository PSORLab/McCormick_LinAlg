module Check_BLAS_MC

    using Compat
    using Compat.Test
    using StaticArrays
    using EAGO #Should be dependent in LinAlg So idk why MC not defined
    using McCormick_LinAlg
    #import McCormick_LinAlg.XSCAL #shouldnt be necessary with proper export

    #using...

#    include("dotproduct.jl")
    include("xscal.jl")
end
