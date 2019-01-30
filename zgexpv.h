#pragma once
#include <iostream>
#include "SLISC/slisc.h"
#include "SLISC/sparse.h"
#include <mkl.h>
#include "zgpadm.h"
#include "znchbv.h"
using std::cout; using std::endl;
using namespace slisc;

void ZGEXPV(Int_I n, Int_I m, Doub_I t, const Comp *v, Comp *w, Doub tol, Doub_I anorm,
	Comp *wsp, Int_I lwsp, Int *iwsp, Int_I liwsp, McooComp_I matvec, Int_I itrace, Int_O iflag )
{
	const Int mxstep = 500, mxreject = 0, ideg = 6;
	Doub delta = 1.2, gamma = 1.9;


	const Comp zero = 0., one = 1.;

	Int i, j, k1, mh, mx, iv, ih, j1v, ns, ifree, lfree, iexph,
		ireject, ibrkflag, mbrkdwn, nmult, nreject, nexph, nscale,
		nstep;
	Doub sgn, t_out, tbrkdwn, step_min, step_max, err_loc,
		s_error, x_error, t_now, t_new, t_step, t_old,
		xm, beta, break_tol, p1, p2, p3, eps, rndoff,
		vnorm, avnorm, hj1j, hump, sqr1;
	Comp hij;
	iflag = 0;
	if (lwsp < n*(m + 2) + 5 * SQR(m + 2) + ideg + 1)
		iflag = -1;
	if (liwsp < m + 2)
		iflag = -2;
	if (m >= n || m <= 0)
		iflag = -3;
	if (iflag <= 0)
		error("bad sizes (in input of ZGEXPV)");

	k1 = 2;
	mh = m + 2;
	iv = 1;
	ih = iv + n*(m + 1) + n;
	ifree = ih + mh*mh;
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
	rndoff = eps*anorm;

	break_tol = 1.0e-7;


	sgn = SIGN(1., t);
	cblas_zcopy(n, v, 1, w, 1);
	beta = cblas_dznrm2(n, w, 1);
	vnorm = beta;
	hump = beta;

	sqr1 = sqrt(0.1);
	xm = 1. / Doub(m);
	p2 = tol*pow((m + 1) / 2.72, m + 1)*sqrt(2.*3.14*(m + 1));
	t_new = (1. / anorm)*pow(p2 / (4.*beta*anorm), xm);
	p1 = pow(10., round(log10(t_new) - sqr1) - 1);
	t_new = trunc(t_new / p1 + 0.55) * p1;

	do {
		if (t_now >= t_out)
			break;

		nstep = nstep + 1;
		t_step = MIN(t_out - t_now, t_new);
		p1 = 1. / beta;
		for (i = 1; i <= n; ++i)
			wsp[iv + i - 1] = p1*w[i];

		for (i = 1; i <= mh*mh; ++i)
			wsp[ih + i - 1] = zero;

		j1v = iv + n;
		for (j = 1; j <= m; ++j) {
			nmult = nmult + 1;
			matvec.mul(wsp + j1v - n, wsp + j1v);
			for (i = 1; i <= j; ++i) {
				cblas_zdotc_sub(n, wsp + iv + (i - 1)*n, 1, wsp + j1v, 1, &hij);
				Comp temp = -hij;
				cblas_zaxpy(n, &temp, wsp + iv + (i - 1)*n, 1, wsp + j1v, 1);
				wsp[ih + (j - 1)*mh + i - 1] = hij;
			}
			hj1j = cblas_dznrm2(n, wsp + j1v, 1);

			if (hj1j <= break_tol) {
				cout << "happy breakdown: mbrkdwn = " << j << "h = " << hj1j << endl;
				k1 = 0;
				ibrkflag = 1;
				mbrkdwn = j;
				tbrkdwn = t_now;
				t_step = t_out - t_now;
				goto 300;
			}
			wsp[ih + (j - 1)*mh + j] = Comp(hj1j);
			cblas_zdscal(n, 1. / hj1j, wsp + j1v, 1);
			j1v = j1v + n;
		}

		nmult = nmult + 1;
		matvec.mul(wsp + j1v - n, wsp + j1v);
		avnorm = cblas_dznrm2(n, wsp + j1v, 1);

		/*300*/
		wsp[ih + m*mh + m + 1] = one;

		ireject = 0;
		/*401*/
		while (true) {
			nexph = nexph + 1;
			mx = mbrkdwn + k1;
			if (ideg != 0) {
				ZGPADM(ideg, mx, sgn*t_step, wsp + ih, mh,
					wsp + ifree, lfree, iwsp, iexph, ns, iflag);
				iexph = ifree + iexph - 1;
				nscale = nscale + ns;
			}
			else {
				iexph = ifree;
				for (i = 1; i <= mx; ++i)
					wsp[iexph + i - 1] = zero;
				wsp[iexph] = one;
				ZNCHBV(mx, sgn*t_step, wsp + ih, mh, wsp + iexph, wsp + ifree + mx);
			}
			/*402*/

			if (k1 == 0) {
				err_loc = tol;
			}
			else {
				p1 = abs(wsp[iexph + m])   * beta;
				p2 = abs(wsp[iexph + m + 1]) * beta * avnorm;
				if (p1 > 10.*p2) {
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
					cout << "t_step = " << t_old << endl;
					cout << "err_loc = " << err_loc << endl;
					cout << "err_required =" << delta*t_old*tol << endl;
					cout << "stepsize rejected, stepping down to:" << t_step << endl;
				}
				ireject = ireject + 1;
				nreject = nreject + 1;
				if (mxreject != 0 && ireject > mxreject) {
					cout << "Failure in ZGEXPV: ---" << endl;
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
		hij = Comp(beta);
		cblas_zgemv(CblasColMajor, CblasNoTrans, n, mx, &hij, wsp + iv, n, wsp + iexph, 1, &zero, w, 1);
		beta = cblas_dznrm2(n, w, 1);
		hump = MAX(hump, beta);

		t_new = gamma * t_step * pow(t_step*tol / err_loc, xm);
		p1 = pow(10., round(log10(t_new) - sqr1) - 1);
		t_new = trunc(t_new / p1 + 0.55) * p1;

		err_loc = MAX(err_loc, rndoff);

		t_now = t_now + t_step;

		if (itrace != 0) {
			cout << "integration " << nstep << "---------------------------------" << endl;
			cout << "scale-square = " << ns << endl;
			cout << "step_size = " << t_step << endl;
			cout << "err_loc   = " << err_loc << endl;
			cout << "next_step = " << t_new << endl;
		}

		step_min = MIN(step_min, t_step);
		step_max = MAX(step_max, t_step);
		s_error = s_error + err_loc;
		x_error = MAX(x_error, err_loc);

	} while (mxstep == 0 || nstep < mxstep);

		iflag = 1;

		/*500*/

		iwsp[1] = nmult;
		iwsp[2] = nexph;
		iwsp[3] = nscale;
		iwsp[4] = nstep;
		iwsp[5] = nreject;
		iwsp[6] = ibrkflag;
		iwsp[7] = mbrkdwn;

		wsp[1] = Comp(step_min);
		wsp[2] = Comp(step_max);
		wsp[3] = Comp(0.);
		wsp[4] = Comp(0.);
		wsp[5] = Comp(x_error);
		wsp[6] = Comp(s_error);
		wsp[7] = Comp(tbrkdwn);
		wsp[8] = Comp(sgn*t_now);
		wsp[9] = Comp(hump / vnorm);
		wsp[10] = Comp(beta / vnorm);
	}