	subroutine xmodelvel(r,val)
c  This program returns PREM velocity for a given depth
c  by the input coefficients of cubic polynomials
c  Input:  	r     ---  radius
c  Output:
c		val   ---  Initial PREM model values
c			   1. ph  2. pv  3. sh  4. sv  5. eta
c
c   Output:	modelvel --- velocity
        implicit double precision (a-h,o-z)
        common/coeff/coef(4,8,20)
	common/path$/delt(2,50000,2),nfin(2),kdep !la 1600->50000
	common/layr/nl(20),xb(20),xt(20),ifanis,nplay
	common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
	common/ray$/pmin,pmax,nps(20,2),npsrc(2),nn,iup,jref,npleg,nsleg
	common/anisopath/theta_path(50000),rad_path(50000),ntheta !la 3200->50000
        common/kernel/sum(362,20,10),ifsplit,ifdiff
        common/isot/isotrp
	dimension val(5)
c
      do j=1, numlyr
	      if(r.ge.xb(j).and.r.lt.xt(j)) then
		      iq = j
		      y=r/rnorm
		      ierror = 0
		      goto 999
	      endif
      enddo
999   if(ierror.eq.1) then
	      write(*,*) 'error in the raypath in getgreen.f'
	      stop
      endif
      val(2)=coef(1,2,iq)+y*(coef(2,2,iq)+y*(coef(3,2,iq)+y*coef(4,2,iq)))
      val(4)=coef(1,3,iq)+y*(coef(2,3,iq)+y*(coef(3,3,iq)+y*coef(4,3,iq)))
      if(.not.isotrp) then
         val(1)=coef(1,6,iq)+y*(coef(2,6,iq)+y*(coef(3,6,iq)+y*coef(4,6,iq)))
         val(3)=coef(1,7,iq)+y*(coef(2,7,iq)+y*(coef(3,7,iq)+y*coef(4,7,iq)))
         val(5)=coef(1,8,iq)+y*(coef(2,8,iq)+y*(coef(3,8,iq)+y*coef(4,8,iq)))
      else
        val(1)=val(2)
        val(3)=val(4)
        val(5)=1.d0
      end if
      return
      end
