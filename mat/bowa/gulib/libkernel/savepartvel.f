      subroutine savepartvel
c
c  saves the partial derivatives and velocities along a path into 
c  the appropriate arrays in a common block.
c
      implicit double precision(a-h,o-z)
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      common/parm/rho,vpv,vsv,vph,vsh,eta,s(5),r
      common/parm1/xa,xc,xn,xl,xf,tdif(5)
      common/isot/isotrp
      common/partvel/pathder(5,50000),pathvel(5,50000)
      common/anisopath/theta_path(50000),rad_path(50000),ntheta
      common/layr/nl(20),xb(20),xt(20),ifanis,nplay
c      dimension pathder(5,3200), pathvel(5,3200) ! commented out by lapo 22.1.2010

      do i=1, ntheta
	r0=rad_path(i)
        do ilay=1, numlyr
	   if(r0.gt.xb(ilay).and.r0.le.xt(ilay)) then
		n0 = ilay
		goto 33
	   endif
      	enddo
33 	call qtau(r0,n0,q0)
	if(q0.ge.0.0) then
c... if SV or SH reflect at fluid boundary, then q2=-1.d0
      		qq0=dsqrt(q0)
      		call tder(r0,qq0)
      		do ipar=1, 5
			pathder(ipar,i)=tdif(ipar)
      		enddo
	else
      		do ipar=1, 5
			pathder(ipar,i)=0.d0
      		enddo				
	endif
	pathvel(1,i)=vph
	pathvel(2,i)=vpv
	pathvel(3,i)=vsh
	pathvel(4,i)=vsv
	pathvel(5,i)=rho
      enddo
      return
      end
