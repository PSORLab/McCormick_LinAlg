@reexport module BLAS_MC
#Not being called right now. Difficulty implementing multiple modules
using StaticArrays, LinearAlgebra

import Base: +, -, *, /, convert, in, isempty, one, zero, real, eps, max, min,
             abs, inv, exp, exp2, exp10, expm1, log, log2, log10, log1p, acosh,
             sqrt, sin, cos, tan, min, max, sec, csc, cot, ^, step, sign, intersect

import IntervalArithmetic: dist, mid, pow, +, -, *, /, convert, in, isempty,
                           one, zero, real, eps, max, min, abs, exp,
                           expm1, log, log2, log10, log1p, sqrt, ^,
                           sin, cos, tan, min, max, sec, csc, cot, step,
                           sign, dist, mid, pow, Interval, interval, sinh, cosh,
                           ∩, IntervalBox, pi_interval, bisect, isdisjoint, length
#export forward operators

export MC, cc, cv, Intv, lo, hi,  cc_grad, cv_grad, cnst, +, -, *, /, convert,
       one, zero, dist, real, eps, mid, exp, exp2, exp10, expm1, log, log2,
       log10, log1p, acosh, sqrt, sin, cos, tan, min, max, sec, csc, cot, ^,
       abs, step, sign, pow, in, isempty, intersect, length

export multtemp, sqrtemp, XSCAL, DOT

export seed_gradient, IntervalType, set_mc_differentiability!, set_multivar_refine!,
       set_tolerance!, set_iterations!

#include all McCormick Arithmatic operations
include("DefaultAttributes\\Constants.jl")
include("DefaultAttributes\\inner_utilities.jl")
include("DefaultAttributes\\Type.jl")
include("DefaultAttributes\\inner_utilities.jl")
include("DefaultAttributes\\Powers.jl")
include("DefaultAttributes\\MultiplicationTemp.jl")
include("DefaultAttributes\\inner_utilities.jl")
#include all BLAS_MC implementations

include("Lin_Functions\\DOT.jl")
include("Lin_Functions\\XSCAL.jl")

#No LAPACK functions yet

end
