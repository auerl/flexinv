!=======================================================================
!
! From rfischer@seismology.harvard.edu Wed Jul 23 15:30:22 1997
! Retrieved by Bob Fischer from: http://www.netlib.org/linalg/lsqr
! From arpa!sol-michael.stanford.edu!mike 5 May 89 23:53:00 PDT
!
      SUBROUTINE LSQR  ( M, N, DAMP,&
                         LENIW, LENRW, IW, rw,&
                         U, V, W, X, SE,&
                         ATOL, BTOL, CONLIM, ITNLIM, NOUT,&
         ISTOP, ITN, ANORM, ACOND, RNORM, ARNORM, XNORM,&
             indx,values,mpoin,nonz,nelrhs)

      integer*8 nonz,mpoin(0:m)
      dimension values(nonz)
      integer*8, dimension(:) :: indx(nonz)


      INTEGER            M, N, LENIW, LENRW, ITNLIM, NOUT, ISTOP, ITN
      INTEGER            IW(LENIW)

      real*4  RW(LENRW), U(M), V(N), W(N), X(N), SE(N),&
                         ATOL, BTOL, CONLIM, DAMP,&
                         ANORM, ACOND, RNORM, ARNORM, XNORM

!-----------------------------------------------------------------------
!     Intrinsics and local variables

      INTRINSIC          ABS, MOD, SQRT
      INTEGER            I, NCONV, NSTOP

      real*4 dnrm2

      real*4   ALFA, BBNORM, BETA, BNORM,&
                         CS, CS1, CS2, CTOL, DAMPSQ, DDNORM, DELTA,&
                         GAMMA, GAMBAR, PHI, PHIBAR, PSI,&
                         RES1, RES2, RHO, RHOBAR, RHBAR1, RHBAR2,&
                         RHS, RTOL, SN, SN1, SN2,&
                         T, TAU, TEST1, TEST2, TEST3,&
                         THETA, T1, T2, T3, XXNORM, Z, ZBAR

      PARAMETER        ( ZERO = 0.,  ONE = 1. )

      CHARACTER*16       ENTER, EXIT
      CHARACTER*60       MSG(0:7)

      DATA               ENTER /' Enter LSQR.    '/,&
                         EXIT  /' Exit  LSQR.    '/

      DATA               MSG &
       / 'The exact solution is  X = 0',&
         'Ax - b is small enough, given ATOL, BTOL',&
         'The least-squares solution is good enough, given ATOL',&
         'The estimate of cond(Abar) has exceeded CONLIM',&
         'Ax - b is small enough for this machine',&
         'The least-squares solution is good enough for this machine',&
         'Cond(Abar) seems to be too large for this machine',&
         'The iteration limit has been reached' /
!-----------------------------------------------------------------------


!     Initialize.

      IF (NOUT .GT. 0)&
         WRITE(NOUT, 1000) ENTER, M, N, DAMP, ATOL, CONLIM, BTOL, ITNLIM
      ITN    =   0
      ISTOP  =   0
      NSTOP  =   0
      CTOL   =   ZERO
      IF (CONLIM .GT. ZERO) CTOL = ONE / CONLIM
      ANORM  =   ZERO
      ACOND  =   ZERO
      BBNORM =   ZERO
      DAMPSQ =   DAMP**2
      DDNORM =   ZERO
      RES2   =   ZERO
      XNORM  =   ZERO
      XXNORM =   ZERO
      CS2    = - ONE
      SN2    =   ZERO
      Z      =   ZERO

      DO 10  I = 1, N
         V(I)  =  ZERO
         X(I)  =  ZERO
        SE(I)  =  ZERO
   10 CONTINUE

!     Set up the first vectors U and V for the bidiagonalization.
!     These satisfy  BETA*U = b,  ALFA*V = A(transpose)*U.

      ALFA   =   ZERO
      BETA   =   DNRM2 ( M, U, 1 )

      IF (BETA .GT. ZERO) THEN
         CALL DSCAL ( M, (ONE / BETA), U, 1 )
         CALL APROD ( 2, M, N, V, U, LENIW, LENRW, IW, RW ,&
             indx,values,mpoin,nonz)
         ALFA   =   DNRM2 ( N, V, 1 )
      END IF

      IF (ALFA .GT. ZERO) THEN
         CALL DSCAL ( N, (ONE / ALFA), V, 1 )
         CALL DCOPY ( N, V, 1, W, 1 )
      END IF

      ARNORM =   ALFA * BETA
      IF (ARNORM .EQ. ZERO) GO TO 800

      RHOBAR =   ALFA
      PHIBAR =   BETA
      BNORM  =   BETA
      RNORM  =   BETA

      IF (NOUT   .GT.  0  ) THEN
         IF (DAMPSQ .EQ. ZERO) THEN
             WRITE(NOUT, 1200)
         ELSE
             WRITE(NOUT, 1300)
         END IF
         TEST1  = ONE
         TEST2  = ALFA / BETA
         WRITE(NOUT, 1500) ITN, X(1), RNORM, TEST1, TEST2
         WRITE(NOUT, 1600)
      END IF

!     ------------------------------------------------------------------
!     Main iteration loop.
!     ------------------------------------------------------------------
  100 ITN    = ITN + 1
	print*,'Iteration:',itn
        print*,'Rnorm/Anorm:',RNORM,ANORM

!     Perform the next step of the bidiagonalization to obtain the
!     next  BETA, U, ALFA, V.  These satisfy the relations
!                BETA*U  =  A*V  -  ALFA*U,
!                ALFA*V  =  A(transpose)*U  -  BETA*V.

      CALL DSCAL ( M, (- ALFA), U, 1 )
!--------verified that calling aprod with m or nelrhs changes nothing. LB 2.2012
      CALL APROD ( 1, M, N, V, U, LENIW, LENRW, IW, RW,&
             indx,values,mpoin,nonz)

      BETA   =   DNRM2 ( M, U, 1 )
      BBNORM =   BBNORM  +  ALFA**2  +  BETA**2  +  DAMPSQ

      IF (BETA .GT. ZERO) THEN
         CALL DSCAL ( M, (ONE / BETA), U, 1 )
         CALL DSCAL ( N, (- BETA), V, 1 )
         CALL APROD ( 2, M, N, V, U, LENIW, LENRW, IW, RW ,&
             indx,values,mpoin,nonz)
         ALFA   =   DNRM2 ( N, V, 1 )
         IF (ALFA .GT. ZERO) THEN
            CALL DSCAL ( N, (ONE / ALFA), V, 1 )
         END IF
      END IF

!     Use a plane rotation to eliminate the damping parameter.
!     This alters the diagonal (RHOBAR) of the lower-bidiagonal matrix.

      RHBAR2 = RHOBAR**2  +  DAMPSQ
      RHBAR1 = SQRT( RHBAR2 )
      CS1    = RHOBAR / RHBAR1
      SN1    = DAMP   / RHBAR1
      PSI    = SN1 * PHIBAR
      PHIBAR = CS1 * PHIBAR

!     Use a plane rotation to eliminate the subdiagonal element (BETA)
!     of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.

      RHO    =   SQRT( RHBAR2  +  BETA**2 )
      CS     =   RHBAR1 / RHO
      SN     =   BETA   / RHO
      THETA  =   SN * ALFA
      RHOBAR = - CS * ALFA
      PHI    =   CS * PHIBAR
      PHIBAR =   SN * PHIBAR
      TAU    =   SN * PHI

!     Update  X, W  and the standard error estimates.

      T1     =   PHI   / RHO
      T2     = - THETA / RHO
      T3     =   ONE   / RHO

      DO 200  I =  1, N
         T      =  W(I)
         X(I)   =  T1*T  +  X(I)
         W(I)   =  T2*T  +  V(I)
         T      = (T3*T)**2
         SE(I)  =  T     +  SE(I)
         DDNORM =  T     +  DDNORM
  200 CONTINUE

!     Use a plane rotation on the right to eliminate the
!     super-diagonal element (THETA) of the upper-bidiagonal matrix.
!     Then use the result to estimate  norm(X).

      DELTA  =   SN2 * RHO
      GAMBAR = - CS2 * RHO
      RHS    =   PHI    - DELTA * Z
      ZBAR   =   RHS    / GAMBAR
      XNORM  =   SQRT( XXNORM    + ZBAR **2 )
      GAMMA  =   SQRT( GAMBAR**2 + THETA**2 )
      CS2    =   GAMBAR / GAMMA
      SN2    =   THETA  / GAMMA
      Z      =   RHS    / GAMMA
      XXNORM =   XXNORM + Z**2

!     Test for convergence.
!     First, estimate the norm and condition of the matrix  Abar,
!     and the norms of  rbar  and  Abar(transpose)*rbar.

      ANORM  =   SQRT( BBNORM )
      ACOND  =   ANORM * SQRT( DDNORM )
      RES1   =   PHIBAR**2
      RES2   =   RES2  +  PSI**2
      RNORM  =   SQRT( RES1 + RES2 )
      ARNORM =   ALFA  * ABS( TAU )

!     Now use these norms to estimate certain other quantities,
!     some of which will be small near a solution.

      TEST1  =   RNORM /  BNORM
      TEST2  =   ZERO
      IF (RNORM .GT. ZERO) TEST2 = ARNORM / (ANORM * RNORM)
      TEST3  =   ONE   /  ACOND
      T1     =   TEST1 / (ONE  +  ANORM * XNORM / BNORM)
      RTOL   =   BTOL  +  ATOL *  ANORM * XNORM / BNORM

!     The following tests guard against extremely small values of
!     ATOL, BTOL  or  CTOL.  (The user may have set any or all of
!     the parameters  ATOL, BTOL, CONLIM  to zero.)
!     The effect is equivalent to the normal tests using
!     ATOL = RELPR,  BTOL = RELPR,  CONLIM = 1/RELPR.

      T3     =   ONE + TEST3
      T2     =   ONE + TEST2
      T1     =   ONE + T1
      IF (ITN .GE. ITNLIM) ISTOP = 7
      IF (T3  .LE. ONE   ) ISTOP = 6
      IF (T2  .LE. ONE   ) ISTOP = 5
      IF (T1  .LE. ONE   ) ISTOP = 4

!     Allow for tolerances set by the user.

      IF (TEST3 .LE. CTOL) ISTOP = 3
      IF (TEST2 .LE. ATOL) ISTOP = 2
      IF (TEST1 .LE. RTOL) ISTOP = 1
!     ==================================================================

!     See if it is time to print something.

      IF (NOUT  .LE.  0       ) GO TO 600
      IF (N     .LE. 40       ) GO TO 400
      IF (ITN   .LE. 10       ) GO TO 400
      IF (ITN   .GE. ITNLIM-10) GO TO 400
      IF (MOD(ITN,10) .EQ. 0  ) GO TO 400
      IF (TEST3 .LE.  2.0*CTOL) GO TO 400
      IF (TEST2 .LE. 10.0*ATOL) GO TO 400
      IF (TEST1 .LE. 10.0*RTOL) GO TO 400
      IF (ISTOP .NE.  0       ) GO TO 400
      GO TO 600

!     Print a line for this iteration.

  400 WRITE(NOUT, 1500) ITN, X(1), RNORM, TEST1, TEST2, ANORM, ACOND
      IF (MOD(ITN,10) .EQ. 0) WRITE(NOUT, 1600)
!     ==================================================================

!     Stop if appropriate.
!     The convergence criteria are required to be met on  NCONV
!     consecutive iterations, where  NCONV  is set below.
!     Suggested value:  NCONV = 1, 2  or  3.

  600 IF (ISTOP .EQ. 0) NSTOP = 0
      IF (ISTOP .EQ. 0) GO TO 100
      NCONV  =   1
      NSTOP  =   NSTOP + 1
      IF (NSTOP .LT. NCONV  .AND.  ITN .LT. ITNLIM) ISTOP = 0
      IF (ISTOP .EQ. 0) GO TO 100
!     ------------------------------------------------------------------
!     End of iteration loop.
!     ------------------------------------------------------------------


!     Finish off the standard error estimates.

      T    =   ONE
      IF (M      .GT.   N )  T = M - N
      IF (DAMPSQ .GT. ZERO)  T = M
      T    =   RNORM / SQRT( T )

      DO 700  I = 1, N
         SE(I)  = T * SQRT( SE(I) )
  700 CONTINUE

!     Print the stopping condition.

  800 IF (NOUT .GT. 0) THEN
         WRITE(NOUT, 2000) EXIT, ISTOP, ITN,&
                           EXIT, ANORM, ACOND,&
                           EXIT, RNORM, ARNORM,&
                           EXIT, BNORM, XNORM
         WRITE(NOUT, 3000) EXIT, MSG(ISTOP)
      END IF

  900 RETURN

!     ------------------------------------------------------------------
 1000 FORMAT(// 1P, A, '  Least-squares solution of  A*x = b'&
          / ' The matrix  A  has', I7, ' rows   and', I7, ' columns'&
          / ' The damping parameter is         DAMP   =', E10.2&
          / ' ATOL   =', E10.2, 15X,        'CONLIM =', E10.2&
          / ' BTOL   =', E10.2, 15X,        'ITNLIM =', I10)
 1200 FORMAT(// '   Itn       x(1)           Function',&
         '     Compatible   LS        Norm A    Cond A' /)
 1300 FORMAT(// '   Itn       x(1)           Function',&
         '     Compatible   LS     Norm Abar Cond Abar' /)
 1500 FORMAT(1P, I6, 2E17.9, 4E10.2)
 1600 FORMAT(1X)
 2000 FORMAT(/ 1P, A, 6X, 'ISTOP =', I3,   16X, 'ITN    =', I9&
             /     A, 6X, 'ANORM =', E13.5, 6X, 'ACOND  =', E13.5&
             /     A, 6X, 'RNORM =', E13.5, 6X, 'ARNORM =', E13.5,&
             /     A, 6X, 'BNORM =', E13.5, 6X, 'XNORM  =', E13.5)
 3000 FORMAT( A, 6X, A )
!     ------------------------------------------------------------------
!     End of LSQR
      END SUBROUTINE LSQR
!=======================================================================




!=======================================================================
      subroutine  dcopy(n,dx,incx,dy,incy)
!
!     copies a vector, x, to a vector, y.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!      double precision dx(*),dy(*)
      real dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
!
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      return
      end subroutine dcopy
!=======================================================================



!=======================================================================
      subroutine  dscal(n,da,dx,incx)
!
!     scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!      double precision da,dx(*)
      real*4 da,dx(*)
      integer i,incx,m,mp1,n,nincx
!
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end subroutine dscal
!=======================================================================



!=======================================================================
!      DOUBLE PRECISION FUNCTION DNRM2 ( N, X, INCX )
      real FUNCTION DNRM2 ( N, X, INCX )
!     .. Scalar Arguments ..
      INTEGER                           INCX, N
!     .. Array Arguments ..
!      DOUBLE PRECISION                  X( * )
      real*4               X( * )
!     ..
!
!  DNRM2 returns the euclidean norm of a vector via the function
!  name, so that
!
!     DNRM2 := sqrt( x'*x )
!
!
!
!  -- This version written on 25-October-1982.
!     Modified on 14-October-1993 to inline the call to DLASSQ.
!     Sven Hammarling, Nag Ltd.
!
!
!     .. Parameters ..
!      DOUBLE PRECISION      ONE         , ZERO
      PARAMETER           ( ONE = 1., ZERO = 0. )
!     .. Local Scalars ..
      INTEGER               IX
!      DOUBLE PRECISION      ABSXI, NORM, SCALE, SSQ
      real*4     ABSXI, NORM, SCALE, SSQ
!     .. Intrinsic Functions ..
      INTRINSIC             ABS, SQRT
!     ..
!     .. Executable Statements ..
      IF( N.LT.1 .OR. INCX.LT.1 )THEN
         NORM  = ZERO
      ELSE IF( N.EQ.1 )THEN
         NORM  = ABS( X( 1 ) )
      ELSE
         SCALE = ZERO
         SSQ   = ONE
!        The following loop is equivalent to this call to the LAPACK
!        auxiliary routine:
!        CALL DLASSQ( N, X, INCX, SCALE, SSQ )
!
         DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
            IF( X( IX ).NE.ZERO )THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI )THEN
                  SSQ   = ONE   + SSQ*( SCALE/ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SSQ   = SSQ   +     ( ABSXI/SCALE )**2
               END IF
            END IF
   10    CONTINUE
         NORM  = SCALE * SQRT( SSQ )
      END IF
!
      DNRM2 = NORM
      RETURN
!
!     End of DNRM2.
!
      END FUNCTION DNRM2
!=======================================================================



!=======================================================================
	subroutine aprod(mode,m,n,x,y,lenin,lenva,iw,aty,&
             indx,values,mpoin,nonz)
        integer*8 nonz,mpoin(0:m),j
	integer iw(lenin)
	real*4 x(n),y(m),aty(n)
	dimension values(nonz)
	integer*8, dimension(:) :: indx(nonz)
	if(mode.eq.1)then
! compute y=y+A*x
	do k=1,m
	   pro=0.
	   do j=mpoin(k-1)+1,mpoin(k)
	      pro=pro+values(j)*x(indx(j))
	   enddo
	   y(k)=y(k)+pro
	enddo
	return
	elseif(mode.eq.2)then
	do i=1,n
	   aty(i)=0.
	enddo
! compute x=x+(A^t)*y
	do k=1,m
	   do j=mpoin(k-1)+1,mpoin(k)
	      aty(indx(j))=aty(indx(j))+values(j)*y(k)
	   enddo
	enddo
	do i=1,n
	   x(i)=x(i)+aty(i)
	enddo
	return
	else
	print*,'error: mode=',mode
	stop
	endif
	end subroutine aprod
!=======================================================================



!=======================================================================
	function isqre(lat,lon,nsqrs,nsqtot,nlatzones,n,eq_incr)
!
!       finds the index of the square where (lat,lon) is
!
	real*4 lat,lon,loc_incr
	dimension nsqrs(nlatzones),nsqtot(nlatzones+1)
	lazone=(90.-lat)/eq_incr+1
	if((90.-lat).gt.180.)lazone=nlatzones
	if((90.-lat).gt.181.)stop "problems in function isqre"
	if(lazone.gt.nlatzones)then
	   print*,"problems in function isqre, latitude",lazone,lat
	   stop
	endif
	if(lon.lt.0.)lon=360.+lon
	if(lon.eq.360.)lon=0.
	loc_incr=360./float(nsqrs(lazone))
	isqre=(lon/loc_incr)+1
	isqre=isqre+nsqtot(lazone)
	if(isqre.gt.n)then
	   print*,"problems in function isqre, longitude",isqre,n,lon,loc_incr,lazone
	   stop
	endif
	return
	end function isqre
!=======================================================================




!=======================================================================
	subroutine damp_norm(weight,nnn,nelp,values,indx,&
       mpoin,t,m,nmantle,nonz,n1)
!
! Norm damping in the regular grid case
!        
        integer*8 nonz,nnn,mpoin(0:m)
	dimension t(m)
	dimension values(nonz)
	integer*8, dimension(:) :: indx(nonz)
	do k=1,N1
	   nnn=nnn+1
	   values(nnn)=WEIGHT
	   indx(nnn)=k+nmantle
	   nnn=nnn+1  ! STH LOOKS WEIRD! BUG? NEED TO CHECK! NOT REALLY USING IT ANYWAYS
	   nelp=nelp+1
	   mpoin(nelp)=nnn
	   t(nelp)=0.
	enddo
	return
	end subroutine damp_norm
!=======================================================================



!=======================================================================
	subroutine damp_norm_adapt(weightvec,nnn,nelp,values,indx,&
       mpoin,t,m,nmantle,nonz,N1)
!
! Norm damping in the variable grid case
! Slightly changed by Ludwig 11/2012
!        
        integer*8 nonz,nnn,mpoin(0:m)
        real*4 weightvec(1000000) ! just set size apriori
	dimension t(m)
	dimension values(nonz)
	integer*8, dimension(:) :: indx(nonz)

	do k=1,N1
	   nnn=nnn+1
	   values(nnn)=weightvec(k)
           print*,weightvec(k)
	   indx(nnn)=k+nmantle
	   nnn=nnn+1 ! THIS LOOKS LIKE A BUG, NOT USING THIS SUBROUTINE ANYHOW, NEED TO CHECK
	   nelp=nelp+1
	   mpoin(nelp)=nnn
	   t(nelp)=0.
	enddo
	
        return
	end subroutine damp_norm_adapt
!=======================================================================




!=======================================================================
	subroutine damp_difference
!
!       Enforces minimal difference between the two inversion parameters
!       Routine is based on Lapos "dampaniso" and was rewritten by Ludwig
!       in May 2014 to allow for a depth variable parameter difference 
!       damping. Works for variable and regular grids.
!
        use flexinv_module

        integer k,l,totad
        
        ie=nelm ! set local ie to global nelm
        ir=nrhs ! set local row index ir to global one icol

        ! In the case of a regular grid parameterization
        if (iadapt.eq.0) then
           do l=1,nlay
              do k=((l-1)*n0+1),(l*n0)
                 !--------------- vph-vpv
                 ie=ie+1
                 val(ie)=wddampvec(l,1)
                 ind(ie)=k
                 ie=ie+1
                 val(ie)=-1.*wddampvec(l,1)
                 ind(ie)=k+n1
                 ir=ir+1
                 pnt(ir)=ie
                 rhs(ir)=0.d0
                 if (npar.eq.4) then
                    !--------------- vsh-vsv
                    ie=ie+1
                    val(ie)=wddampvec(l,2)
                    ind(ie)=k+2*n1
                    ie=ie+1
                    val(ie)=-1.*wddampvec(l,2)
                    ind(ie)=k+3*n1
                    ir=ir+1
                    pnt(ir)=ie
                    rhs(ir)=0.d0
                 end if
              enddo
           enddo
        ! In the case of a variable grid parameterization
        else if (iadapt.eq.1) then
           totad=0
           do l=1,nlay
              do k=totad+1,totad+adpx1layer(l)
                 !--------------- vph-vpv
                 ie=ie+1
                 val(ie)=wddampvec(l,1)
                 ind(ie)=k
                 ie=ie+1
                 val(ie)=-1.*wddampvec(l,1)
                 ind(ie)=k+n1
                 ir=ir+1
                 pnt(ir)=ie
                 rhs(ir)=0.d0
                 if (npar.eq.4) then
                    !--------------- vsh-vsv
                    ie=ie+1
                    val(ie)=wddampvec(l,2)
                    ind(ie)=k+2*n1
                    ie=ie+1
                    val(ie)=-1.*wddampvec(l,2)
                    ind(ie)=k+3*n1
                    ir=ir+1
                    pnt(ir)=ie
                    rhs(ir)=0.d0
                 end if
              end do
              totad=totad+adpx1layer(l)
           enddo
        end if
        
        ! restore global indices
        nelm   = ie
        nrhs   = ir

	end subroutine damp_difference
!=======================================================================


!=======================================================================
	subroutine prescribe_aniso(preani,preani_model,nnn,nelp,&
        values,indx,mpoin,rhs,m,n_1,NONZ)
!
!  Forces the anisotropic part of the solution in the case of a xi-voigt_average
!  paramterization to be equal to a input hypothesis
!
        integer*8 nonz,nnn,mpoin(0:m)
	dimension rhs(m)
	dimension values(nonz)
	integer*8, dimension(:) :: indx(nonz)
        real*4 xival,anival,dmy
        character*140 preani_model
        open(747,file=trim(preani_model))
        print*,"Prescribing anisotropy"
	do k=1+n_1,n_1*2     
           print*,"Parameter: ",k-n_1
           read(747,*) dmy,xival
           print*,anival
           anival=xival
	   nnn=nnn+1
	   values(nnn)=preani
	   indx(nnn)=k
	   nelp=nelp+1
	   mpoin(nelp)=nnn
	   rhs(nelp)=anival*preani
	enddo
        close(747)
	return
	end subroutine prescribe_aniso
!=======================================================================



!=======================================================================
	function superisqre(lat,lon,nsqrs,nsqtot,nlatzones,n,eq_incr,nlatzomax)
!
!  Finds the index of square where lat, lon is
!
	dimension nsqrs(nlatzomax),nsqtot(nlatzomax+1)
	incr=eq_incr*100.
	lazone=(9000-lat)/incr+1
	if(lazone.gt.nlatzones)lazone=nlatzones
	llon=lon
	if(llon.lt.0)llon=36000+llon
        superisqre=(llon*nsqrs(lazone))/36000+1
	superisqre=superisqre+nsqtot(lazone)
	if(superisqre.gt.n)superisqre=n
	return
	end function superisqre
!=======================================================================



!=======================================================================
	SUBROUTINE RANGE(NSQ,XLAMIN,XLAMAX,XLAMID,XLOMIN,XLOMAX,&
      	XLOMID,nsqrs,nsqtot,nlatzones,n,eq_incr)
! 
! FINDS THE COORDINATE RANGE OF SQUARE NUMBER 'NSQ'
!
        PARAMETER(NLATZOMAX=180.)
	DIMENSION NSQRS(NLATZONES),NSQTOT(NLATZONES)
	LAZONE=2
	DO WHILE (NSQ.GT.NSQTOT(LAZONE))
	   LAZONE=LAZONE+1
	ENDDO
	LAZONE=LAZONE-1
	NNSQ=NSQ-NSQTOT(LAZONE)
	XLAMIN=90.-LAZONE*eq_incr
	XLAMAX=XLAMIN+eq_incr
	XLAMID=XLAMIN+(eq_incr*0.5)
	GRSIZE=360./NSQRS(LAZONE)
	XLOMAX=NNSQ*GRSIZE
	XLOMIN=XLOMAX-GRSIZE
	XLOMID=XLOMAX-(GRSIZE*0.5)
	RETURN
	END SUBROUTINE RANGE! end of subroutine RANGE
!=======================================================================




!=======================================================================
	subroutine coordsuper(nlatzomax,nbloc,blocla,bloclo,nsqrs,nlatzones,eq_incr)
!
!  given a cell index on the Earth's surface, finds longitude and latitude
!  (NOT colatitude) of its center.
!
	dimension nsqrs(nlatzomax)

	ntot=0
! loop(s) over all the blocks
	do 500 ila=1,nlatzones
! increment latitude
	   rlati=90.-(eq_incr*(ila-1))
! calculate increment in longitude for this band
	   rinlo=(360./nsqrs(ila))
	   do 400 isq=1,nsqrs(ila)
	      rlong=(360./nsqrs(ila))*(isq-1)
	      ntot=ntot+1
	      if(ntot.eq.nbloc)then
	         bloclo=rlong+(rinlo/2.)
	         blocla=rlati-(eq_incr/2.)
	         goto 600
	      endif
400	   continue
500	continue
600	return
	end subroutine coordsuper
!=======================================================================



!=======================================================================
        subroutine contribution_ata(row,index,ata,d,b,n,nata,nz)
!
!  given a row of a and corresponding datum augments ata accordingly
!  j it's unimportant for ata which row, just index of value in row important  
!
        real*8 ata(nata),b(n)
        real row(n),d
        integer*8 index(n)
        do i=1,nz
           do j=i,nz              
              if((index(j).eq.0).or.(index(i).eq.0))then
                 print*,"ERROR",index(j),index(i),i,j
                 stop
              end if
              if(index(j).ge.index(i))then
                 ind=(((index(j)-1)*index(j))/2)+index(i)
              else
                 ind=(((index(i)-1)*index(i))/2)+index(j)
              endif
              ata(ind)=ata(ind)+row(i)*row(j)
              b(index(i))=b(index(i))+(row(i)*d)
           enddo
           ! add to rhs vector the contribution of j-th row:
        enddo
!        print*,"-----"
        return
        end subroutine contribution_ata
!=======================================================================





!=======================================================================
        subroutine contribution_ata2(row,index,ata,d,b,n,nata,nz)
!
!  given a row of a and corresponding datum augments ata accordingly
!  j it's unimportant for ata which row, just index of value in row important  
!
        real*8 ata(nata),b(n)
        real row(n),d
        integer*8 index(n)
        do i=1,nz
!           print*,row(i),d,nz,index(i)
           do j=i,nz              
              if((index(j).eq.0).or.(index(i).eq.0))then
                 print*,"ERROR",index(j),index(i),i,j
                 stop
              end if
              if(index(j).ge.index(i))then
                 ind=(((index(j)-1)*index(j))/2)+index(i)
              else
                 ind=(((index(i)-1)*index(i))/2)+index(j)
              endif
              ata(ind)=ata(ind)+row(i)*row(j)
              b(index(i))=b(index(i))+(row(i)*d)
           enddo
        enddo
        return
        end subroutine contribution_ata2
!=======================================================================






!=======================================================================
      subroutine damp_roughness 

!
! Written by Ludwig in Winter 2012, efficient method to set up the
! "finite difference" roughness damping operator. Replaces Julias
! Spherical Harmonic based approach
!

                                                               
      use flexinv_module


      ! Local variables      
      integer, parameter :: nperad=30
      integer, dimension(nperad) :: ib1,ib2,ib3
      integer, dimension(nperad) :: ib4,ib5,ib6
      integer, dimension(nlaym) :: layers

      integer, allocatable, dimension(:,:) :: level !(n0max,nlaym)
      real, dimension(4) :: weight

      real, allocatable, dimension(:,:) :: xlamin,xlamax
      real, allocatable, dimension(:,:) :: xlomin,xlomax

      integer i,j,k,l ! i layer index
      integer u,v,p,z ! u npar index
      integer ir,zo,ia,ie ! ie=nelm index, ir=nrhs index
      integer kk1,kk2,kk3
      integer kk4,kk5,kk6
      integer nadpar,clay

      real lomin,lomax,lamin,lamax
      real romin,romax,ramin,ramax
      real lama,lami,loma,lomi
      real dmy,weightvv,weightg
      real xx1,xx2,xx3,xx4
      real bloclo,blocla

      character*3 i3      
      character*1 dmy3

      ! Allocate memory
      allocate(xlamin(n0max,nlaym),xlamax(n0max,nlaym))
      allocate(xlomin(n0max,nlaym),xlomax(n0max,nlaym))
      allocate(level(n0max,nlaym),valb(ntemp),indb(ntemp))

      ! read in grid information
      nadpar=0
      if (iadapt.eq.1) then
         open(10,file=trim(gridinfo))
         open(20,file=trim(gridinfo)//".sh")
         open(30,file=trim(gridinfo)//".lay")
         open(40,file=trim(gridinfo)//".gmt")
         do i=1,nlay
            read(30,*),layers(i)
            nadpar=nadpar+layers(i)
            write(i3,"(i3.3)")i         
            ! separate grid in layers and export them to output folder
            open(50,file='./grid.layer.'//i3//'.gmt')
            open(60,file='./grid.layer.'//i3//'.sh')
            do j=1,layers(i)
               read(10,*)  ,dmy,level(j,i)
               read(20,*)  ,dmy,xlamin(j,i),xlamax(j,i),xlomin(j,i),xlomax(j,i)
               write(60,*) ,dmy,xlamin(j,i),xlamax(j,i),xlomin(j,i),xlomax(j,i)
               read(40,*)  ,dmy3
               write(50,*) ,dmy3
               read(40,*)  ,xx1,xx2
               write(50,*) ,xx1,xx2
               read(40,*)  ,xx1,xx2
               write(50,*) ,xx1,xx2
               read(40,*)  ,xx1,xx2
               write(50,*) ,xx1,xx2
               read(40,*)  ,xx1,xx2
               write(50,*) ,xx1,xx2
            end do
            close(50)
            close(60)
         end do
         close(10)
         close(20)
         close(30)
         close(40)
      elseif (iadapt.eq.0) then         
         do i=1,nlay
            do j=1,n0
               call range(j,lami,lama,blocla,lomi,loma,bloclo,&
                         nsqrs,nsqtot,nlatzones,n0,eqincr)
               level(j,i)  = 1
               layers(i)   = n0
               xlamax(j,i) = lama
               xlamin(j,i) = lami
               xlomax(j,i) = loma
               xlomin(j,i) = lomi
            end do
            nadpar=nadpar+n0
         end do
      end if
      close(99)

!
! major loop over npar 1 and 2, number of layers and 
! number of pixels in each layer
!

      ie=nelm ! set local ie to global nelm
      ir=nrhs ! set local row index ir to global one icol

      if (iata.eq.1) then
         ie=0
         ir=0
      end if

      do u=1,npar

         ! individual dmaping factors for each param
         clay=0
         do i=1,nlay

            ! hardcoded, smaller weights for smaller voxels
            weight(1) = wgradh(i,u)
            weight(2) = wgradh(i,u)*0.9
            weight(3) = wgradh(i,u)*0.8 
            weightvv  = wgradv(i)

            print*,"Damping layer",i,"with:",weight(1),weightvv
            do j=1,layers(i)

               ! identify the current pixel index
               ie=ie+1
               ia=(u-1)*nadpar+clay+j
               
               if (ia.gt.(nadpar*npar)) then                 
                  stop "ERROR: Undefined block index! Stopping ..."
               end if
               
               ind(ie)=ia
               val(ie)=0.
               
               ! second loop, over all
               ramin=xlamin(j,i)
               ramax=xlamax(j,i)
               romin=xlomin(j,i)
               romax=xlomax(j,i)
               
               kk1=0
               kk2=0
               kk3=0
               kk4=0
               kk5=0
               kk6=0
               
               ! 2nd loop over all pixels in matrix
               do z=1,3
                  if(z.eq.1)then
                     zo=i ! this is the curr
                  elseif(z.eq.2)then
                     if(i.eq.1)then
                        goto 999 ! skip this round of top layer, as we are in first layer
                     else
                        zo=i-1
                     endif
                  elseif(z.eq.3)then
                     if(i.eq.nlay)then
                        goto 999 ! skip this round of top layer, as we are in first layer
                     else
                        zo=i+1
                     endif
                  endif
                  
                  do k=1,layers(zo)
                     
                     ! where are we right now
                     lamin=xlamin(k,zo)
                     lamax=xlamax(k,zo)
                     lomin=xlomin(k,zo)
                     lomax=xlomax(k,zo)
                     
                     if(z.eq.3)then
                        if((lamin.ge.ramin).and.(lamax.le.ramax).and.&
                                  (lomin.ge.romin).and.(lomax.le.romax))then
                           kk5=kk5+1
                           ib5(kk5)=(u-1)*nadpar+(clay+layers(i))+k
                        elseif((ramin.ge.lamin).and.(ramax.le.lamax).and.&
                                  (romin.ge.lomin).and.(romax.le.lomax))then
                           kk5=kk5+1
                           ib5(kk5)=(u-1)*nadpar+(clay+layers(i))+k
                        endif
                     elseif(z.eq.2)then
                        if((lamin.ge.ramin).and.(lamax.le.ramax).and.&
                                  (lomin.ge.romin).and.(lomax.le.romax))then
                           kk6=kk6+1
                           ib6(kk6)=(u-1)*nadpar+(clay-layers(i-1))+k   ! removed a bug here i -> i-1
                        elseif((ramin.ge.lamin).and.(ramax.le.lamax).and.&
                                  (romin.ge.lomin).and.(romax.le.lomax))then
                           kk6=kk6+1
                           ib6(kk6)=(u-1)*nadpar+(clay-layers(i-1))+k   ! removed a bug here i -> i-1
                        endif
                     elseif(z.eq.1)then
                        if(lamin.eq.ramax)then ! south adjacent
                           if((lomin.ge.romin).and.(lomax.le.romax))then
                              kk1=kk1+1
                              ib1(kk1)=(u-1)*nadpar+clay+k
                           elseif((romin.ge.lomin).and.(romax.le.lomax))then
                              kk1=kk1+1
                              ib1(kk1)=(u-1)*nadpar+clay+k
                           elseif((lomin.lt.romax).and.(lomax.gt.romax))then
                              kk1=kk1+1
                              ib1(kk1)=(u-1)*nadpar+clay+k
                           elseif((lomin.lt.romin).and.(lomax.gt.romin))then
                              kk1=kk1+1
                              ib1(kk1)=(u-1)*nadpar+clay+k
                           endif
                        elseif(lamax.eq.ramin)then ! north adjacent
                           if((lomin.ge.romin).and.(lomax.le.romax))then
                              kk2=kk2+1
                              ib2(kk2)=(u-1)*nadpar+clay+k
                           elseif((romin.ge.lomin).and.(romax.le.lomax))then
                              kk2=kk2+1
                              ib2(kk2)=(u-1)*nadpar+clay+k
                           elseif((lomin.lt.romax).and.(lomax.gt.romax))then
                              kk2=kk2+1
                              ib2(kk2)=(u-1)*nadpar+clay+k
                           elseif((lomin.lt.romin).and.(lomax.gt.romin))then
                              kk2=kk2+1
                              ib2(kk2)=(u-1)*nadpar+clay+k
                           endif
                        elseif(lomin.eq.romax)then ! east adjacent
                           if((lamin.ge.ramin).and.(lamax.le.ramax))then
                              kk3=kk3+1
                              ib3(kk3)=(u-1)*nadpar+clay+k
                           elseif((ramin.ge.lamin).and.(ramax.le.lamax))then
                              kk3=kk3+1
                              ib3(kk3)=(u-1)*nadpar+clay+k
                           endif
                        elseif(lomax.eq.romin)then ! west adjacent
                           if((lamin.ge.ramin).and.(lamax.le.ramax))then
                              kk4=kk4+1
                              ib4(kk4)=(u-1)*nadpar+clay+k
                           elseif((ramin.ge.lamin).and.(ramax.le.lamax))then
                              kk4=kk4+1
                              ib4(kk4)=(u-1)*nadpar+clay+k
                           endif
                           ! special treatment at boundary
                        elseif((lomax.eq.360.).and.(romin.eq.0.))then ! west
                           if((lamin.ge.ramin).and.(lamax.le.ramax))then
                              kk4=kk4+1
                              ib4(kk4)=(u-1)*nadpar+clay+k
                           elseif((ramin.ge.lamin).and.(ramax.le.lamax))then
                              kk4=kk4+1
                              ib4(kk4)=(u-1)*nadpar+clay+k
                           endif
                        elseif((lomin.eq.0.).and.(romax.eq.360.))then ! east
                           if((lamin.ge.ramin).and.(lamax.le.ramax))then
                              kk3=kk3+1
                              ib3(kk3)=(u-1)*nadpar+clay+k
                           elseif((ramin.ge.lamin).and.(ramax.le.lamax))then
                              kk3=kk3+1
                              ib3(kk3)=(u-1)*nadpar+clay+k
                           endif
                        endif
                     endif
                  end do ! end of loop over pixels in same layer
999               continue ! gets here in case of top or bottom layer
               end do ! loop over top mid and bot layer
               
               ! Increment the matrix
               do p=1,kk1
                  val(ie)=val(ie)+real((1./real(kk1))*weight(level(j,i)))            
               end do
               do p=1,kk2
                  val(ie)=val(ie)+real((1./real(kk2))*weight(level(j,i)))                        
               end do
               do p=1,kk3
                  val(ie)=val(ie)+real((1./real(kk3))*weight(level(j,i)))            
               end do
               do p=1,kk4
                  val(ie)=val(ie)+real((1./real(kk4))*weight(level(j,i)))            
               end do
               do p=1,kk5
                  val(ie)=val(ie)+real((1./real(kk5))*weight(level(j,i))*weightvv)            
               end do
               do p=1,kk6
                  val(ie)=val(ie)+real((1./real(kk6))*weight(level(j,i))*weightvv)            
               end do
 
               ! Increment the matrix
               do v=1,kk1
                  ie=ie+1
                  val(ie)=real((-1./real(kk1))*weight(level(j,i)))
                  ind(ie)=ib1(v)
                  if(ib1(v).gt.(nadpar*npar))then
                     print*,"undefined index"
                     stop
                  endif
               end do
               do v=1,kk2
                  ie=ie+1
                  val(ie)=real((-1./real(kk2))*weight(level(j,i)))
                  ind(ie)=ib2(v)
                  if(ib2(v).gt.(nadpar*npar))then
                     print*,"undefined index"
                     stop
                  endif
               end do
               do v=1,kk3
                  ie=ie+1
                  val(ie)=real((-1./real(kk3))*weight(level(j,i)))
                  ind(ie)=ib3(v)
                  if(ib3(v).gt.(nadpar*npar))then
                     print*,"undefined index"
                     stop
                  endif
               end do
               do v=1,kk4
                  ie=ie+1
                  val(ie)=real((-1./real(kk4))*weight(level(j,i)))
                  ind(ie)=ib4(v)
                  if(ib4(v).gt.(nadpar*npar))then
                     print*,"undefined index"
                     stop
                  endif
               end do
               do v=1,kk5
                  ie=ie+1
                  val(ie)=real((-1./real(kk5))*weight(level(j,i))*weightvv)
                  ind(ie)=ib5(v)
                  if(ib5(v).gt.(nadpar*npar))then
                     print*,"undefined index"
                     stop
                  endif
               end do
               do v=1,kk6
                  ie=ie+1
                  val(ie)=real((-1./real(kk6))*weight(level(j,i))*weightvv)
                  ind(ie)=ib6(v)
                  if(ib6(v).gt.(nadpar*npar))then
                     print*,"undefined index"
                     stop
                  endif
               end do              
               
               ! increment t and pnt vectors
               ir=ir+1
               rhs(ir)=0.
               pnt(ir)=ie
               
            end do
            clay = clay + layers(i) ! we go to next layer
         end do
      end do

! NOTE ATA/CHOLESKY MODE IS DEFUNCTIONAL
!      if (iata.eq.1) then
!         do u=1,n
!            if (mod(u,int(n/100)).eq.0) then
!               print*,"Added,",int(100*u/n)+1, "% of damping rows to ata"
!           end if
!           irec=0
!            do p=pnt(u-1)+1,pnt(u)
!               irec=irec+1
!               valb(irec)=val(p)
!               indb(irec)=ind(p)
!            end do
!            print*,irec,valb(1:irec),indb(1:irec)
!            call contribution_ata2(valb(1:irec),indb(1:irec),&
!                 ata,rhs(u),atd,n,natamax,irec)
!         end do
!      end if

      deallocate(xlamin,xlamax)
      deallocate(xlomin,xlomax)
      deallocate(level,valb,indb)

      ! Set back global index to local index
      nelm   = ie
      nrhs   = ir
      end subroutine damp_roughness
!=======================================================================







!=======================================================================
        subroutine param(eqincr,nsqrs,nsqtot,nlatzones,numto,iswit,&
                         refgrid,nlatzomax,iadapt)
!
! this version is parameterized exactly as julias version
! and allows to distinguish between adaptive and regular
! parameterization, L.A. Aug. 2012
!
        integer nlatzones0
        integer nsqrs(nlatzomax),nsqtot(nlatzomax+1)
        integer nsqrs0(nlatzomax),nsqtot0(nlatzomax+1)
	parameter(pi=3.1415926536)


        ! this is ludwigs preferred regular parameterization
        if(iadapt.eq.0)then
           print*,"Working with Ludwigs prefered regular grid parameterization"
           numto=0
           colat=-eqincr/2.
           do k=1,nlatzones
              colat=colat+eqincr
              theta=(colat/180.)*pi
              ! for this latitudinal zone, compute number of blocks (nsqrs)
              deltalon=eqincr/(sin(theta))
              nsqrs(k)=(360./deltalon)+1
              ! if requested, correct nsqrs(k) so the grid is compatible to reference grid
              if(iswit.eq.1)then
                 if(360./nsqrs(k).ge.refgrid)then
 100             if((mod(360./nsqrs(k),refgrid).ne.0).or.&
                    (mod(nsqrs(k),2).ne.0))then
                    nsqrs(k)=nsqrs(k)+1
                    goto 100
                 else
                 endif
              elseif(360./nsqrs(k).lt.refgrid)then
 101             if((mod(refgrid,360./nsqrs(k)).ne.0).or.&
                    (mod(nsqrs(k),2).ne.0))then
                    nsqrs(k)=nsqrs(k)-1
                    goto 101
                 else
                 endif
              endif
           else ! modified by lud
              if(mod(nsqrs(k),2).ne.0)nsqrs(k)=nsqrs(k)-1
	   endif
         if(mod(nsqrs(k),2).ne.0)then
              stop "nsqrs has to be even"
         endif
         ! print regular grid parameterization
         print*,k,nsqrs(k)
         nsqtot(k)=numto
         numto=numto+nsqrs(k)
         enddo
	 nsqtot(nlatzones+1)=numto         


         ! this is julias adaptive grid, which is slighl different from mine
         elseif(iadapt.eq.1)then
           print*,"Working with julias prefered adaptive grid parameterization"
            n1layer0=0
            nlatzones0=180./refgrid
            colat=-refgrid/2.

            do k=1,nlatzones0
               colat=colat+refgrid
               theta=(colat/180.)*pi
               ! for this latitudinal zone, compute number of blocks (nsqrs)
               deltalon=refgrid/(sin(theta))
               nsqrs0(k)=(360./deltalon)+1
               if(mod(nsqrs0(k),2).ne.0)nsqrs0(k)=nsqrs0(k)-1
               ! if requested, correct nsqrs(k) so the grid is compatible to reference grid
               if(iswit.eq.1)then
                  if(360./nsqrs0(k).ge.refgrid)then
 102              if(mod(360./nsqrs0(k),refgrid).ne.0)then
                     nsqrs0(k)=nsqrs0(k)+1
                     goto 102
                  else
                  endif
               elseif(360./nsqrs0(k).lt.refgrid)then
 103              if(mod(refgrid,360./nsqrs0(k)).ne.0)then
                     nsqrs0(k)=nsqrs0(k)-1
                     goto 103
                  else
                  endif
               endif
	       endif
               
               if(mod(nsqrs0(k),2).ne.0)stop "nsqrs has to be even"
               ! print rough parameterization acc. julia
               print*,k,nsqrs0(k)
               nsqtot0(k)=n1layer0
               n1layer0=n1layer0+nsqrs0(k)
            enddo

            ! define fine grid
            fact=refgrid/eqincr
            n1layer=0
            print*,"ratio of coarser to finest grid is",fact
            do k=1,nlatzones
               k0=((k-1)/fact)+1
               nsqrs(k)=nsqrs0(k0)*fact
               nsqtot(k)=n1layer
               n1layer=n1layer+nsqrs(k)
               print*, "k=",k,"k0=",k0," nsqrs=", nsqrs(k),nsqrs0(k0)," lon=", 360./nsqrs(k)
               print*, "latitudinal zone =",k, "lon=", 360./nsqrs(k), "#pixels=", nsqrs(k) 
            enddo
            nsqtot(nlatzones+1)=n1layer
            print*,'Number of pixel with finest parameterization: n1layer=',n1layer
            numto=n1layer
            
         endif
       
         return
         end subroutine param
!=======================================================================

