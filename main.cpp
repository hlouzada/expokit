#include <iostream>
#include "SLISC/slisc.h"
#include "SLISC/sparse.h"
#include "SLISC/disp.h"
#include "SLISC/arithmatic.h"
#include <mkl.h>
#include "zhexpv.h"

using namespace slisc;
using std::cout; using std::endl;

void test_ZHEXPV()
{
	Int Nsize = 6;
	McooComp A;
	A.resize(100);
	A.resize(Nsize, Nsize);
	A.set(1, 0, 0);
	A.set(2, 1, 1);
	A.set(3, 2, 2);
	A.set(4, 3, 3);
	A.set(5, 4, 4);
	A.set(6, 5, 5);
	A.set(Comp(0., 1), 0, 5);
	A.set(Comp(0., -1.), 5, 0);
	VecComp v(Nsize), w(Nsize);
	linspace(v, 1, Nsize);
	w = 0.;
	Doub t = 1;
	Doub tol = 1e-10;
	Doub anorm = 1.;
	Int N_krylov = 3;
	Int lwsp = SQR(Nsize*(N_krylov + 2) + 5 * (N_krylov + 2)) + 7;
	VecComp wsp(lwsp);
	Int liwsp = N_krylov + 2;
	VecInt iwsp(liwsp);
	Int itrace = 1;
	Int iflag;
	ZHEXPV(Nsize, N_krylov, t, v.ptr(), w.ptr(), tol, anorm,
		wsp.ptr(), lwsp, iwsp.ptr(), liwsp, A, itrace, iflag);
	disp(w);
}

int main()
{
	test_ZHEXPV();

	/*Int Nsize = 2;
	McooComp A;
	A.resize(100);
	A.resize(Nsize, Nsize);
	A.set(1, 0, 0);
	A.set(2, 1, 1);
	A.set(Comp(0., 1.), 1, 0);
	A.set(Comp(0., -1.), 0, 1);
	VecComp v(Nsize), w(Nsize);
	linspace(v, 1, 2);
	w = 0.;
	mul(w.ptr(), A, v.ptr());
	disp(w);*/
}
