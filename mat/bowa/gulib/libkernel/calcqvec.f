      subroutine calcqvec(r0,q)
c
c  Compute dq/dv as given by Woodhouse 1981.
c  Input:	r0=radius
c  Output:	q(1..5)=derivative dq/dv, etc.
c
      implicit double precision(a-h,o-z)
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      common/parm/rho,vpv,vsv,vph,vsh,eta,s(5),r
      common/parm1/xa,xc,xn,xl,xf,tdif(5)
      common/isot/isotrp
      common/anisopath/theta_path(50000),rad_path(50000),ntheta !la 3200->50000
      common/layr/nl(20),xb(20),xt(20),ifanis,nplay
      dimension q(1)


c *** Identify in which layer we are right now
      do ilay=1, numlyr
	   if(r0.gt.xb(ilay).and.r0.le.xt(ilay)) then
		n0 = ilay
		goto 33
	   endif
      enddo


c *** Compute the radicand of qtau
33    call qtau(r0,n0,q0)

      if(q0.ge.0.0) then
c *** Note: if SV or SH reflect at fluid boundary, then q2=-1.d0

	    qq0=dsqrt(q0)  ! Note: qtau only gives the radicand

c *** Compute (analytically) first derivatives of qq0 wrt to model parameters
	    call tder(r0,qq0)
          
c *** Assign the output q to tdif (which comes from tder)
          do ipar=1, 5
		q(ipar)=tdif(ipar)
          enddo
      else
          do ipar=1, 5
		q(ipar)=0.d0
          enddo				
      endif

      return
      end
