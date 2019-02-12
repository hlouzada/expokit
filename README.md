# Expokit (C++ translation) 
[Expokit](https://www.maths.uq.edu.au/expokit/download.html) is a numerical library that calculates the matrix exponential or it's operation on a vector.

This project is a manual translation of some expokit subroutines to C++. With a trivial modification, one can change them to C as well.

This project uses some utilities in the [SLISC](https://github.com/MacroUniverse/SLISC) library, mainly for custom type definition and a COO Sparse matrix class that supports matrix-vector multiplication.

## Translation Details
* The original code it translated as directly as possible except the following changes.
* Since C++ arrays are 0-based, some variables always have value that is -1 than it's fortran counterpart. Note that this include some function arguments. See the comments before each function for a complete list.
* since `goto` statement is deprecated, it's replaced with other structures.
* BLAS and LAPACK subroutines are replaced with corresponding CBLAS and LAPACKE from Intel's MKL library (currently free). With minor modifications, other CBLAS and LAPACKE may be used.

## Tests
* Hermitian matrices up to 40-by-40 are tested with ZGEXPV and ZHEXPV, with results 2-3 digits less than machine precision.
