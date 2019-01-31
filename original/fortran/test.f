*
*---  sample program illustrating the use of ZGEXPV and ZHEXPV ...
*     Hermitian problem (Example 6.2 in the Expokit report) ...
*
      implicit none
      external zgcoov

*---  matrix data ...
*---  BEWARE: these values must match those in zgmatv.f
      integer n, nz, nmax, nzmax
      parameter( nmax = 5500, nzmax = 50000 )
      integer ia(nzmax), ja(nzmax)
      complex*16 a(nzmax)
      common /CMAT/ a, ia, ja, nz, n

*---  arguments variables ...
      integer m, mmax, lwsp, liwsp
      parameter( mmax = 50 )
      parameter( lwsp = nmax*(mmax+2)+5*(mmax+2)**2+7, liwsp = mmax+2 )
      integer iwsp(liwsp)
      double precision t, tol, anorm, s1, s2
      complex*16 v(nmax), w(nmax), wsp(lwsp)

      integer i, j, ii, kk, nnz, itrace, iflag
      complex*16 ZERO, ONE
      parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )

      intrinsic ABS, CMPLX, CONJG, DBLE
*
*---  Executable statements ...

*---  load a symmetric pattern ...
      n = 20
      nz = 3*n - 2;

*---  for the purpose of the experiments, expand to COOrdinates ...
      kk = 0
      do ii = 1,n
         kk = kk + 1
         ia(kk) = ii; ja(kk) = ii; a(kk) = 1d0;
         v(ii) = ii
      enddo
      do ii = 1,n-1
         kk = kk + 1
         ia(kk) = ii; ja(kk) = ii + 1; a(kk) = (0.1d0);
      enddo
      do ii = 1,n-1
         kk = kk + 1
         ia(kk) = ii + 1; ja(kk) = ii; a(kk) = (0.1d0);
      enddo

      do i = 1,n
         wsp(i) = ZERO
      enddo
      do i = 1,nz
         wsp(ia(i)) = wsp(ia(i)) + ABS( a(i) )
      enddo
      anorm = wsp(1)
      do i = 2,n
         if ( anorm.lt.DBLE(wsp(i)) ) anorm =  wsp(i)
      enddo

      t = 1.0d0
      tol = 0d0
      m = 10
      itrace = 1

      call ZGEXPV( n, m, t,v,w, tol, anorm,
     .             wsp,lwsp, iwsp,liwsp, zgcoov, itrace, iflag )
      print*, 'done testing!'
      end
