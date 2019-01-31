#include <iostream>
#include "SLISC/slisc.h"
#include "SLISC/sparse.h"
#include "SLISC/disp.h"
#include "SLISC/arithmatic.h"
#include <mkl.h>
#include "zhexpv.h"
#include "zgexpv.h"

using namespace slisc;
using std::cout; using std::endl;

Doub norm_inf(McooComp_I A)
{
	VecDoub abs_sum(A.nrows(), 0.);
	for (Int i = 1; i < A.nnz(); ++i) {
		abs_sum(A.row(i)) += abs(A(i));
	}
	return max(abs_sum);
}

void test_ZHEXPV()
{
	// === params ===========
	Int n = 20; // matrix size
	Int m = 10; // # krylov basis
	Doub t = 1; // time
	Doub tol = 0; // error tol
	// ======================

	Int i;
	McooComp A;
	A.resize(3*n - 2);
	A.resize(n, n);
	
	for (i = 0; i < n; ++i) {
		A.push(1., i, i);
	}
	for (i = 0; i < n - 1; ++i) {
		A.push(0.1, i, i + 1);
	}
	for (i = 0; i < n - 1; ++i) {
		A.push(0.1, i + 1, i);
	}

	VecComp v(n), w(n, 0.);
	linspace(v, 1, n);
	
	Doub anorm = norm_inf(A);

	Int lwsp = MAX(10, SQR(n*(m + 2) + 5 * (m + 2)) + 7);
	VecComp wsp(lwsp);
	Int liwsp = MAX(m + 1, 7);
	liwsp = liwsp + 1;
	VecInt iwsp(liwsp);
	Int itrace = 1;
	Int iflag;
	ZGEXPV(n, m, t, v.ptr(), w.ptr(), tol, anorm,
		wsp.ptr(), lwsp, iwsp.ptr(), liwsp, A, itrace, iflag);
	disp(w, 16);
}

void disp(const Comp *ptr, Int_I n)
{
	Int i;
	for (i = 0; i < 100; ++i) {
		cout << "\n";
	}
	for (i = 0; i < n; ++i) {
		cout << ptr[i] << endl;
	}
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
