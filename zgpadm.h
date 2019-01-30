#pragma once
#include "SLISC/slisc.h"
#include <iostream>
#include <mkl.h>

using namespace slisc;
using std::cout; using std::endl;

// translate as directly as possible for now:
// * any variable values should not change
// * change all variables to lower case
// * don't change BLAC/LAPACK routines
// * all pointer indexing (including +, -) is subtracted by 1 in place

void ZGPADM(Int_I ideg, Int_I m, Doub_I t, const Comp *H, Int_I ldh, Comp *wsp, Int_I lwsp,
	Int *ipiv, Int_O iexph, Int_O ns, Int_O iflag)
{
	Int i, j, k, icoef, mm, ih2, iodd, iused, ifree, iq, ip, iput, iget;
	Doub hnorm;
	Comp cp, cq, scale, scale2, temp;
	const Comp zero = 0., one = 1.;

	mm = m*m;
	iflag = 0;
	if (ldh < m) iflag = -1;
	if (lwsp < 4 * mm + ideg + 1) iflag = -2;
	if (iflag != 0)
		error("bad sizes (in input of ZGPADM)");

	icoef = 1;
	ih2 = icoef + (ideg + 1);
	ip = ih2 + mm;
	iq = ip + mm;
	ifree = iq + mm;

	for (i = 1; i <= m; ++i) {
		wsp[i - 1] = zero;
	}
	for (j = 1; j <= m; ++j) {
		for (i = 1; i <= m; ++i) {
			wsp[i - 1] = wsp[i - 1] + abs(H[i + ldh*(j-1) - 1]);
		}
	}

	hnorm = 0.;
	for (i = 1; i <= m; ++i) {
		hnorm = MAX(hnorm, real(wsp[i - 1])); // not sure what DBLE(wsp(i)) means (wsp(i) is Comp)
	}

	hnorm = abs(t*hnorm);
	if (hnorm == 0.)
		error("Error - null H in input of ZGPADM.");
	ns = MAX(0, Int(log(hnorm) / log(2.)) + 2);
	scale = Comp(t / pow(2, ns), 0.);
	scale2 = scale*scale;

	i = ideg + 1;
	j = 2 * ideg + 1;
	wsp[icoef - 1] = one;
	for (k = 1; k <= ideg; ++k) {
		wsp[icoef + k - 1] = (wsp[icoef + k - 2]*Doub(i - k)) / Doub(k*(j - k));
	}

	cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, m, m, &scale2, H, ldh,
		H, ldh, &zero, wsp + ih2 - 1, m);

	cp = wsp[icoef + ideg - 2];
	cq = wsp[icoef + ideg - 1];
	for (j = 1; j <= m; ++j) {
		for (i = 1; i <= m; ++i) {
			wsp[ip + (j - 1)*m + i - 2] = zero;
			wsp[iq + (j - 1)*m + i - 2] = zero;
		}
		wsp[ip + (j - 1)*(m + 1) - 1] = cp;
		wsp[iq + (j - 1)*(m + 1) - 1] = cq;
	}

	iodd = 1;
	k = ideg - 1;

	/*100*/
	do {
		iused = iodd*iq + (1 - iodd)*ip;
		cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, m, m, &one, wsp + iused - 1, m,
			wsp + ih2 - 1, m, &zero, wsp + ifree - 1, m);
		for (j = 1; j <= m; ++j)
			wsp[ifree + (j - 1)*(m + 1) - 1] = wsp[ifree + (j - 1)*(m + 1) - 1] + wsp[icoef + k - 2];

		ip = (1 - iodd)*ifree + iodd*ip;
		iq = iodd*ifree + (1 - iodd)*iq;
		ifree = iused;
		iodd = 1 - iodd;
		k = k - 1;
	} while (k > 0);


	if (iodd != 0) {
		cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, m, m, &scale, wsp + iq - 1, m,
			H, ldh, &zero, wsp + ifree - 1, m);
		iq = ifree;
	}
	else {
		cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, m, m, &scale, wsp + ip - 1, m,
			H, ldh, &zero, wsp + ifree - 1, m);
		ip = ifree;
	}
	temp = -one;
	cblas_zaxpy(mm, &temp, wsp + ip - 1, 1, wsp + iq - 1, 1);
	iflag = LAPACKE_zgesv(LAPACK_COL_MAJOR, m, m, (MKL_Complex16*)wsp + iq - 1, m, ipiv,
		(MKL_Complex16*)wsp + ip - 1, m);
	if (iflag != 0)
		error("Problem in ZGESV (within ZGPADM)");
	cblas_zdscal(mm, 2., wsp + ip - 1, 1);
	for (j = 1; j <= m; ++j)
		wsp[ip + (j - 1)*(m + 1) - 1] = wsp[ip + (j - 1)*(m + 1) - 1] + one;

	iput = ip;
	if (ns == 0 && iodd != 0) {
		cblas_zdscal(mm, -1., wsp + ip - 1, 1);
	}
	else {
		iodd = 1;
		for (k = 1; k <= ns; ++k) {
			iget = iodd*ip + (1 - iodd)*iq;
			iput = (1 - iodd)*ip + iodd*iq;
			cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, m, m, &one, wsp + iget - 1, m,
				wsp + iget - 1, m, &zero, wsp + iput - 1, m);
			iodd = 1 - iodd;
		}
	}

	/*200*/
	iexph = iput;
}
