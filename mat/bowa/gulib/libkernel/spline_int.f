      subroutine spline_int(eqlat,eqlon,stlat,stlon,bslat,bslon,ngpt,
     +	h,rknt,nknt,numsplit,isplayer,xtu,iray)
c
c  Computes the contribution of a spline function over
c  a ray path.
c  Input:
c	 eqlat ---  Earthquake source lat
c	 eqlon ---  Earthquake source lon
c	 stlat ---  Earthquake receiver lat
c	 stlon ---  Earthquake receiver lon
c	 bslat ---  lat array of horizontal splines
c	 bslon ---  lon array of horizontal splines
c	 ngpt  ---  number of horizontal grid points
c	 h     ---  bucket size (arc length of nonzero B-spline contribution)
c	 rknt  ---  radial splines
c	 nknt  ---  number of radial splines
c	 numsplit  ---  number of split in the radial knots of desired 3-D model
c	 isplayer   ---  buffer containing top and bottom index of split mantle
c	 theta_path   ---  array of running angles along ray path
c	 rad_path   ---  array of radius along raypath
c	 nthetaeta     ---  total number of path segment
c	 xtu     ---  turning radius of the ray
c	 
c  Output:
c 	 arow --- a row of A matrix 
c 	 indarow --- index of non-zero elements of A 
c
      implicit double precision (a-h,o-z)

c *** basic parameters, naming is obsolete
      parameter(maxfreq=8)
      parameter(maxnode=maxfreq*maxfreq*10+2)
      parameter(maxrknot=50)
      parameter(maxparm=maxnode*maxrknot)
      data pi /3.141592653579/
      data rad/57.2957795130824/
      data reprad /57.295779579/    


c... A-matrix 
      common/partvel/pathder(5,3200),pathvel(5,3200)
      common/amatrix/arow(maxparm),indarow(maxparm),nonzero_a
      common/invopt/numinvvar,invparlist(5)
      common/radcontrib/qseg(5,3200),arowk(5,maxrknot)
      common/anisopath/theta_path(3200),rad_path(3200),ntheta
      dimension qrad(3,3200),fwork(3,3200)
      dimension qth(3,3200),qrad_l(3,3200),qth_l(3,3200),xleng(3200)
      dimension bslat(1),bslon(1)
      dimension isplayer(maxrknot,1)
      dimension rknt(maxrknot)
      dimension transmat(9)
c--- save the segment of path that a spline contributes to----
c---
c      dimension rseg(3200),qseg(5,3200),thseg(3200), 
      dimension  xlseg(3200)
c-------------------------------------------------------------

      hpi=pi*0.5d0
      twoh=2.d0*h
      twohsq=twoh*twoh
      call conv2geocen(eqlat, eqlon, eqtheta, eqphi)
      call conv2geocen(stlat, stlon, sttheta, stphi)
      del = dist(eqtheta, eqphi, sttheta, stphi) ! distance between source + receiver
      call calc_fwdtransmat(eqtheta, eqphi, sttheta, stphi, transmat)
cTEST
	open(88,file="ray2d.txt")
	do i=1,ntheta
	write(88,*)theta_path(i)*180./pi,rad_path(i)
	enddo
	close(88)
c
c... compute length of ray dl=sqrt(dr+dth*r), put in array xleng
      xl=0.0
c      xleng(i)=0.d0 !this looked like a typo - lapo 22.2.2010
      xleng(1)=0.d0
      do i=1,ntheta-1
	 dr0=rad_path(i+1)-rad_path(i)
	 dx0=(theta_path(i+1)-theta_path(i))*(rad_path(i+1)+rad_path(i))*0.5
	 dl=sqrt(dx0*dx0+dr0*dr0)
	 xl=xl+dl
	 xleng(i+1)=xl ! length along ray in km
      enddo
      xleng(ntheta+1)=xleng(ntheta) 
c... cubic spline interpolation of radius and theta with respect to xleng
      call drspln(1, ntheta, theta_path, rad_path, qrad, fwork)
      call drspln(1, ntheta, xleng, rad_path, qrad_l, fwork)
      call drspln(1, ntheta, xleng, theta_path, qth, fwork)
      call drspln(1, ntheta, theta_path, xleng, qth_l, fwork)
     
      do i=1, ngpt
        call conv2geocen(bslat(i), bslon(i), sptheta, spphi)
c... transform spline location to source--receiver coordinates
        call transform(sptheta, spphi, transmat, spth1, spph1)
c... check if spline is within 2*h of equator.
	if(abs(spth1-hpi).le.twoh) then
	   aa=hpi-spth1
	   bb=sqrt(twohsq-aa*aa)	! distance from center of bell
	   xmin=spph1-bb 		! boundaries of the bell
	   xmax=spph1+bb 
	   if(xmax.gt.twopi) then
	   	xmax=xmax-twopi
	   	xmin=xmin-twopi
	   endif
	   mode=1
	   if(ifbetween(0.d0,xmin,del).ne.0.and.ifbetween(0.d0,xmax,del).ne.0) then
c... path passes bell
	   	ii0=-1
	   	ii1=-1
	   else if(ifbetween(0.d0,xmin,del).ne.0.and.xmax.gt.del) then
c... path ends in bell
	   	xmax=del
	   	ii0=-1
	   	ii1=ntheta
	   else if(ifbetween(0.d0,xmax,del).ne.0.and.xmin.lt.0) then
c... path begins in bell
	    	xmin=0.d0
	   	ii0=1
	   	ii1=-1
	   else if(ifbetween(xmin,0.d0,xmax).ne.0.and.ifbetween(xmin,del,xmax).ne.0) then
c... path all in bell
	   	xmin=0.d0
	   	xmax=del
	   	ii0=1
	   	ii1=ntheta
	   else
	   	mode=0
	   endif
c
c ... below finds beginning and ending limits and 
c ... integrate radially
c
           if(mode.ne.0) then
c...below saves index of point at xmin, lower limit
	   	if(ii0.lt.0) then
			do j=1,ntheta
			   if(theta_path(j).ge.xmin) goto 10
			enddo
10			ii0=j
	   	endif
c...below saves index of point at xmax, upper limit
	   	if(ii1.lt.0) then
			do j=ntheta,1,-1
			   if(theta_path(j).le.xmax) goto 20
			enddo
20	  		ii1=j
	   	endif

c
c ...The following sets the entire path for the integration.  The idea is that 
c ...the sum of the int(dq/dv1+dq/dv2)=-(travel time) if this works well.
c ...To test that, uncomment the following four lines and activate the code related
c ...to sum(i) in radcontribute_drdl.f.  My experience is that S and ScS works very
c ...well (sum of kernels agrees with travel time to first decimal place) and 
c ...it is a couple of seconds off for the SS waves.
c		ii0=2
c		ii1=ntheta-1
c		xmin=theta_path(ii0)-(theta_path(ii0)-theta_path(ii0-1))*0.5
c		xmax=theta_path(ii1)+(theta_path(ii1+1)-theta_path(ii1))*0.5
c------------------------------------------------------------------------------------
		indseg=0		
	      	if(ii0.ne.1) then
c
c... account for the extra point (xmin) at the beginning. Add points at xmin and 
c... at theta_path(ii0)
c
			indseg=indseg+1
          		xlseg(indseg)=drsple(1,ntheta,theta_path,xleng,qth_l,xmin)
c			indseg=indseg+1
c          		rr=drsple(1,ntheta,theta_path,rad_path,qrad,xmin)
c			dx0=(theta_path(ii0)-xmin)*(rr+rad_path(ii0))*0.5
c			dr0=rad_path(ii0)-rr
c			xlseg(indseg)=sqrt(dx0*dx0+dr0*dr0)
	  	endif
c
c... copy the segment of the ray to the appropriate arrays-----------------------------
c
		do j=0, ii1-ii0
			indseg=indseg+1
			ind0=j+ii0
c			rseg(indseg)=rad_path(ind0)
c			thseg(indseg)=theta_path(ind0)
			xlseg(indseg)=xleng(ind0)
		enddo
c---account for the extra point at the end---xmax---------------------------------------
	   	if(ii1.ne.ntheta) then
			indseg=indseg+1
          		rr=drsple(1,ntheta,theta_path,rad_path,qrad,xmax)
			dx1=(xmax-theta_path(ii1))*(rr+rad_path(ii1))*0.5
			dr1=rad_path(ii1)-rr
			xlseg(indseg)=xleng(ii1)+sqrt(dx1*dx1+dr1*dr1)
	   	endif
		nseg=indseg
c---------------------------------------------------------------------------------------
c  integrate over splines using Gauss-Legendre 5-point integration
c		print*, spth1, spph1, h
c		print*, xleng
	   	call radcontribute_drdl(spth1,spph1,h,xleng,xlseg,nseg,qrad_l,qth,
     +				rknt,nknt,numsplit,isplayer,xtu)
	   	ind1=0
		do ipar=1, numinvvar
			ind=invparlist(ipar)
	   		do ik=1, nknt
				ind2=ind1+i
				if(iray.eq.1) then
					arow(ind2)=arow(ind2)+arowk(ind,ik)
				else
c... subtracting second phase for differential times
					arow(ind2)=arow(ind2)-arowk(ind,ik)
				endif
				ind1=ind1+ngpt
			enddo
	   	enddo
	   endif
	 endif
      enddo
      return
      end
