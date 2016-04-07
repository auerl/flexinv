      subroutine getmod_iso(y,iq,vpeq,vseq)
c
c... computes equivalent velocity from anisotropic velocities
c    
c
      implicit double precision (a-h,o-z)
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      common/layr/nl(20),xb(20),xt(20),ifanis,nplay
      common/isot/isotrp
c
c... get anisotropic velocities
      call getmod(y,iq,rho,vpv,vsv,vph,vsh,eta)
      xc=rho*vpv*vpv	! C (check PREM paper)
      xl=rho*vsv*vsv	! L
      xa=xc		! A = C if isotropic
c
c  below is the isotropic velocity constructed from anisotropic
c  components
c
      if(.not.isotrp) then
        xa=rho*vph*vph	! A
        xn=rho*vsh*vsh	! C
        xf=eta*(xa-2.d0*xl)	! F = eta*(A-2*L)
        xkapa=(4.d0*xa+xc+4.d0*xf-4.d0*xn)/9.d0
        xmu=(xa+xc-2.d0*xf+5.d0*xn+6.d0*xl)/15.d0
        xa=xkapa+4.d0*xmu/3.d0
        xf=xkapa-2.d0*xmu/3.d0
        xn=xmu 
        vph=dsqrt(xa/rho)
        vsh=dsqrt(xn/rho)
	vpeq=vph   ! these are equivalent velocities
	vseq=vsh
      else
	vpeq=vpv
	vseq=vsv
      endif
      return
      end
