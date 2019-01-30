#include <iostream>
#include "SLISC/slisc.h"
#include "SLISC/sparse.h"
#include <mkl.h>
#include "zgpadm.h"
#include "znchbv.h"

using std::cout; using std::endl;
using namespace slisc;

// translate as directly as possible for now:
// * any variable values should not change
// * change all variables to lower case
// * don't change BLAC/LAPACK routines
// * all pointer indexing (including +, -) is subtracted by 1 in place

void ZHEXPV(Int_I n, Int_I m, Doub_I t, const Comp *v, Comp *w, Doub tol, Doub_I anorm,
			Comp *wsp, Int_I lwsp, Int *iwsp, Int_I liwsp, McooComp_I matvec, Int_I itrace, Int_O iflag )
{
	const Int mxstep = 500, mxreject = 0, ideg = 6;
	Doub delta = 1.2, gamma = 0.9;

	const Comp zero = 0., one = 1.;

	Int i, j, k1, mh, mx, iv, ih, j1v, ns, ifree, lfree, iexph,
		ireject, ibrkflag, mbrkdwn, nmult, nreject, nexph, nscale,
		nstep;
	Doub sgn, t_out, tbrkdwn, step_min, step_max, err_loc,
		s_error, x_error, t_now, t_new, t_step, t_old,
		xm, beta, break_tol, p1, p2, p3, eps, rndoff,
		vnorm, avnorm, hj1j, hump, sqr1;
	Comp hjj, temp;

	iflag = 0;
	if ( lwsp < SQR(n*(m+2)+5*(m+2))+ideg+1 ) iflag = -1;
	if (liwsp < m + 2) iflag = -2;
	if (m >= n || m <= 0) iflag = -3;
	if (iflag != 0) error("bad sizes (in input of DHEXPV)");

	k1 = 2;
	mh = m + 2;
	iv = 1;
	ih = iv + n * (m + 1) + n;
	ifree = ih + mh * mh;
	lfree = lwsp - ifree + 1;

	ibrkflag = 0;
	mbrkdwn = m;
	nmult = 0;
	nreject = 0;
	nexph = 0;
	nscale = 0;

	t_out = abs(t);
	tbrkdwn = 0.;
	step_min = t_out;
	step_max = 0.;
	nstep = 0;
	s_error = 0.;
	x_error = 0.;
	t_now = 0.;
	t_new = 0.;

	p1 = 4. / 3.;
	do {
		/*1*/  p2 = p1 - 1.;
		p3 = p2 + p2 + p2;
		eps = abs(p3 - 1.);
	} while (eps == 0.);

	if (tol <= eps) tol = sqrt(eps);
	rndoff = eps * anorm;

	break_tol = 1e-7;

	sgn = t < 0 ? -1 : t > 0 ? 1. : 0.;
	cblas_zcopy(n, v, 1, w, 1);
	beta = cblas_dznrm2(n, w, 1);
	vnorm = beta;
	hump = beta;

	sqr1 = sqrt(0.1);
	xm = 1. / Doub(m);
	p2 = tol * (pow((m + 1) / 2.72, m + 1))*sqrt(2.*3.14*(m + 1));
	t_new = (1. / anorm)*pow((p2 / (4.*beta*anorm)), xm);
	p1 = pow(10., round(log10(t_new) - sqr1) - 1);
	t_new = trunc(t_new / p1 + 0.55) * p1;


	/*100*/ do {
		if (t_now >= t_out)
			break;

		nstep = nstep + 1;
		t_step = MIN(t_out - t_now, t_new);
		beta = cblas_dznrm2(n, w, 1);
		p1 = 1. / beta;
		for (i = 1; i <= n; ++i) {
			wsp[iv + i - 2] = p1 * w[i-1];
		}
		for (i = 1; i <= mh * mh; ++i) {
			wsp[ih + i - 2] = zero;
		}

		j1v = iv + n;
		Bool break_flag = false;
		/*200*/ for (j = 1; j <= m; ++j) {
			nmult = nmult + 1;
			mul(wsp + j1v - n - 1, matvec, wsp + j1v - 1);
			if (j > 1) {
				temp = -wsp[ih + (j - 1)*mh + j - 3];
				cblas_zaxpy(n, &temp, wsp + j1v - 2 * n - 1, 1, wsp + j1v - 1, 1);
			}
				
			cblas_zdotc_sub(n, wsp + j1v - n - 1, 1, wsp + j1v - 1, 1, &hjj);
			temp = -hjj;
			cblas_zaxpy(n, &temp, wsp + j1v - n - 1, 1, wsp + j1v - 1, 1);
			hj1j = cblas_dznrm2(n, wsp + j1v - 1, 1);
			wsp[ih + (j - 1)*(mh + 1) - 1] = hjj;

			if (hj1j <= break_tol) {
				cout << "happy breakdown: mbrkdwn =" << j << " h = " << hj1j << endl;
				k1 = 0;
				ibrkflag = 1;
				mbrkdwn = j;
				tbrkdwn = t_now;
				t_step = t_out - t_now;
				break_flag = true;
				break;
			}
			wsp[ih + (j - 1)*mh + j - 1] = Comp(hj1j);
			wsp[ih + j * mh + j - 2] = Comp(hj1j);
			cblas_zdscal(n, 1. / hj1j, wsp + j1v - 1, 1);
			j1v = j1v + n;
		}

		if (!break_flag) {
			nmult = nmult + 1;
			mul(wsp + j1v - n - 1, matvec, wsp + j1v - 1);
			avnorm = cblas_dznrm2(n, wsp + j1v - 1, 1);
		}

		/*300*/
		wsp[ih + m * mh + m - 2] = zero;
		wsp[ih + m * mh + m] = one;

		ireject = 0;

		while (true) {
			nexph = nexph + 1;
			mx = mbrkdwn + k1;
			if (ideg != 0) {
				ZGPADM(ideg, mx, sgn*t_step, wsp + ih - 1, mh,
					wsp + ifree - 1, lfree, iwsp, iexph, ns, iflag);
				iexph = ifree + iexph - 1;
				nscale = nscale + ns;
			}
			else {
				iexph = ifree;
				for (i = 1; i <= mx; ++i) {
					wsp[iexph + i - 2] = zero;
				}
				wsp[iexph - 1] = one;
				ZNCHBV(mx, sgn*t_step, wsp + ih - 1, mh, wsp + iexph - 1, wsp + ifree + mx - 1);
			}

			/*402*/

			if (k1 == 0) {
				err_loc = tol;
			}
			else {
				p1 = abs(wsp[iexph + m - 1])   * beta;
				p2 = abs(wsp[iexph + m]) * beta * avnorm;
				if (p1 > 10. * p2) {
					err_loc = p2;
					xm = 1. / Doub(m);
				}
				else if (p1 > p2) {
					err_loc = (p1*p2) / (p1 - p2);
					xm = 1. / Doub(m);
				}
				else {
					err_loc = p1;
					xm = 1. / Doub(m - 1);
				}
			}

			if ((k1 != 0) && (err_loc > delta*t_step*tol) &&
				(mxreject == 0 || ireject < mxreject)) {
				t_old = t_step;
				t_step = gamma * t_step * pow(t_step*tol / err_loc, xm);
				p1 = pow(10., round(log10(t_step) - sqr1) - 1);
				t_step = trunc(t_step / p1 + 0.55) * p1;
				if (itrace != 0) {
					cout << "t_step =" << t_old << endl;
					cout << "err_loc =" << err_loc << endl;
					cout << "err_required =" << delta * t_old*tol << endl;
					cout << "stepsize rejected, stepping down to:" << t_step << endl;
				}
				ireject = ireject + 1;
				nreject = nreject + 1;
				if (mxreject != 0 && ireject > mxreject) {
					cout << "Failure in ZHEXPV: ---" << endl;
					cout << "The requested tolerance is too high." << endl;
					cout << "Rerun with a smaller value." << endl;
					iflag = 2;
					return;
				}
				continue;
			}
			break;
		}

		mx = mbrkdwn + MAX(0, k1 - 1);
		hjj = Comp(beta);
		cblas_zgemv(CblasColMajor, CblasNoTrans, n, mx, &hjj,
			wsp + iv - 1, n, wsp + iexph - 1, 1, &zero, w, 1);
		beta = cblas_dznrm2(n, w, 1);
		hump = MAX(hump, beta);

		t_new = gamma * t_step * pow(t_step*tol / err_loc, xm);
		p1 = pow(10., round(log10(t_new) - sqr1) - 1);
		t_new = trunc(t_new / p1 + 0.55) * p1;

		err_loc = MAX(err_loc, rndoff);

		t_now = t_now + t_step;

		if (itrace != 0) {
			cout << "integration" << nstep << "---------------------------------" << endl;
			cout << "scale-square =" << ns << endl;
			cout << "step_size = " << t_step << endl;
			cout << "err_loc   =" << err_loc << endl;
			cout << "next_step =" << t_new << endl;
		}

		step_min = MIN(step_min, t_step);
		step_max = MAX(step_max, t_step);
		s_error = s_error + err_loc;
		x_error = MAX(x_error, err_loc);

	} while (mxstep == 0 || nstep < mxstep);

	iflag = 1;

	/*500*/

	iwsp[0] = nmult;
	iwsp[1] = nexph;
	iwsp[2] = nscale;
	iwsp[3] = nstep;
	iwsp[4] = nreject;
	iwsp[5] = ibrkflag;
	iwsp[6] = mbrkdwn;

	wsp[0] = Comp(step_min);
	wsp[1] = Comp(step_max);
	wsp[2] = Comp(0.);
	wsp[3] = Comp(0.);
	wsp[4] = Comp(x_error);
	wsp[5] = Comp(s_error);
	wsp[6] = Comp(tbrkdwn);
	wsp[7] = Comp(sgn*t_now);
	wsp[8] = Comp(hump / vnorm);
	wsp[9] = Comp(beta / vnorm);
}
