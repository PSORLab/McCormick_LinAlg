# McCormick_LinAlg: Linear Algebra For McCormick Relaxations
McCormick_LinAlg is an auxillary library of linear algebra functions designed for the global optimization solver EAGO (Easy Advanced Global Optimization). For more information about EAGO, click [**here**](https://github.com/PSORLab/EAGO.jl).

## Features
* BLAS styled linear algebra routines meant to handle the McCormick relaxation data type found in EAGO.
* Test files to check correctness and efficiency of each function's implementations
* Estimation methods for the determinant of an interval matrix. **Still in progress**


## Status of Package
 * Majority of methods in the [**BLAS**](http://www.netlib.org/blas/#_optimized_blas_library) library up to TRMV are running.
 * Preconditioned Gaussian Elimination seems to be the most functional and accurate Interval Determinant Estimation method
