	subroutine calc_arclen(iradfun)
c  The following routine calculates the arc length of the
c  ray  and derivative for each segment. 
c  Input:	
c	     iradfun --- choice of radial function
c		kmax --- number of radial B-spline
c  Output:
c		sum  ---  kernel saved in /pathker/ 
c	       rlen  ---  arc length of the ray in each segment
c
	implicit double precision(a-h,o-z)
	integer kmax, iradfun
	real*8 qrad(3, 1600+2), fwork(3,1600+2), b(4)
	dimension qray(4,5,20)
	common/path$/delt(2,800,2),nfin(2),kdep
	common/layr/nl(20),xb(20),xt(20),ifanis,nplay
	common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
	common/ray$/pmin,pmax,nps(20,2),npsrc(2),nn,iup,jref,npleg,nsleg
	common/savepath/theta_path(1600),rad_path(1600)
     +       ,ntheta, nraytp(20),nrayst(20),ileg
        common/invparm/sum(362,20,10), rlen(1600), rker(-1:40),  nspline,
     +	      kmax, krad(1600), rfac(3,1600)
c
c initialize sum()
c	
	if(iradfun.gt.0.and.iradfun.lt.3) then
		do i=1, 5
			do j=1, kmax+1
				do k=1, nspline
					sum(k,j,i) = 0.d0
				enddo
			enddo
		enddo
	endif
c change ray path in terms of radius
	do i=1, ntheta
		rad_path(i) = rad_path(i)*6371.d0
	enddo
c
c Cubic Spline interpolation to find the coeffients of the cubic polynomial
c note that ray path really starts at rad_path(2), not rad_path(1)
c
	rad_path(ntheta+1) = rad_path(ntheta)
	theta_path(ntheta+1) = theta_path(ntheta)
	rr0 = rad_path(2)
	pp0 = theta_path(2)
	do i=2, ntheta
		if(i.ne.(ntheta)) then
			rr1 = rad_path(i+1)
			pp1 = theta_path(i+1)
		else
			rr1 = rr0
			pp1 = pp0
		endif
		krad(i) = getreg(rad
		ireg = krad(i)
		if(ireg.gt.0.and.ireg.le.maxk) then
			call dbspfun(ireg,rad_path(i),rknt,maxk,b)
			vv = 0.5d0/
			rfac(i,1) =vv*b(1) 
			rfac(i,2) =vv*b(2) 
			rfac(i,3) =vv*b(3) 
			rfac(i,4) =vv*b(4) 
		else
			rfac(1,i)=0.d0
			rfac(2,i)=0.d0
			rfac(3,i)=0.d0
			rfac(4,i)=0.d0
		endif
		dr=rr1-rr0
		dphi=(pp1-pp0)*rr1
		rlen(i)=dsqrt(dr*dr+dphi*dphi)
		rr0 = rr1
		pp0 = pp1
	enddo
	return
	end
