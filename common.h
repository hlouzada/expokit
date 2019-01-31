#pragma once
#include <iostream>
#include <mkl.h>

#include "SLISC/arithmatic.h"
#include "SLISC/sparse.h"

// basic types
using slisc::Bool;
using slisc::Int; using slisc::Int_I; using slisc::Int_O;
using slisc::Doub; using slisc::Doub_I;
using slisc::Comp;

// basic inline functions
using slisc::MIN;
using slisc::MAX;
using slisc::SQR;
using slisc::SIGN;

// sparse matrix class
using slisc::McooComp; using slisc::McooComp_I;
using slisc::mul;
