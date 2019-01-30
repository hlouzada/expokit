#include <iostream>
#include "SLISC/slisc.h"
#include "SLISC/disp.h"
#include "SLISC/arithmatic.h"
#include <mkl.h>

using namespace slisc;
using std::cout; using std::endl;

int main()
{
	VecComp x, y;
	linspace(x, Comp(0,1), Comp(19,20), 20);
	disp(x);
	y.resize(20);
	cblas_zcopy(20, x.ptr(), 1, y.ptr(), 1);
	disp(y);
}
