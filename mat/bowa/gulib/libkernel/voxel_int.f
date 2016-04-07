      subroutine voxel_int(eqlat,eqlon,edep,stlat,stlon,bslat,bslon,ngpt,
     +	h,rknt,nknt,numsplit,isplayer,xtu,iray,nsqrs,nsqtot,rbnd,nlatzones,nlay,
     +	nlatzomax,nlaym,eq_incr,n1layer,rowps,indro,irec,ierr,nvx,angle_hitcount)
c
c  Computes the contribution of a voxel function over a ray path.
c  Input:
c	 eqlat ---  Earthquake source lat
c	 eqlon ---  Earthquake source lon
c	 edep  ---  Earthquake depth
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

c *** to output raypaths for diagnostics
      dimension thpath(50000),radpath(50000)

c *** A-matrix 
      common/partvel/pathder(5,50000),pathvel(5,50000) ! la 3200->50000, whole block
      common/amatrix/arow(maxparm),indarow(maxparm),nonzero_a
      common/invopt/numinvvar,invparlist(5)
      common/radcontrib/qseg(5,50000),arowk(5,maxrknot)
      common/anisopath/theta_path(50000),rad_path(50000),ntheta
      dimension qrad(3,50000),fwork(3,50000)
      dimension qth(3,50000),qrad_l(3,50000),qth_l(3,50000),xleng(50000)
      dimension bslat(1),bslon(1)
      dimension isplayer(maxrknot,1)
      dimension rknt(maxrknot)
      dimension transmat(9)

c *** New parameters/arrays for voxel version 
      dimension xlseg(50000)
      parameter(ndim=50000,rearth=6371.d0)
      parameter(nptint=10)    ! default: 10, points inside voxel to evaluate path integral
      parameter(nint=5)       ! default: 3, how many interpolations of ray path
      dimension xray(ndim),yray(ndim),theta_path_deg(ndim)
      dimension xi(ndim*2),yi(ndim*2),delint(ndim*2),zi(ndim*2)
      dimension nsqrs(nlatzomax),nsqtot(nlatzomax+1),rbnd(0:nlaym)
      parameter(toler=0.00005)! default:0.00005 numerical error in pixel coordinates
      dimension raydtmp(ndim),ztmp(ndim)
      dimension rowps(4,50000),indro(50000)
      integer isave(ndim)
      integer nvx
      real*4 angle_hitcount(4,31808*2*50) ! TODO: TOO MANY STATIC ARRAYS

c *** Some parameters
      ierr=0
      hpi=pi*0.5d0
      twoh=2.d0*h
      twohsq=twoh*twoh



c *** Converts from geographical to geocentered coordinates
      call conv2geocen(eqlat, eqlon, eqtheta, eqphi) !
      call conv2geocen(stlat, stlon, sttheta, stphi) ! 
      del = dist(eqtheta, eqphi, sttheta, stphi)     ! dist b/w source + receiver in radians

      
c *** convert epdist in case of major arc 
      if (arc.eq.1) then 
         del=360-del
      endif

      call calc_fwdtransmat(eqtheta, eqphi, sttheta, stphi, transmat) ! do i still need this?


c *** Convert theta_path to degree
	do i=1,ntheta
	   theta_path_deg(i)=theta_path(i)*180.d0/pi
	enddo

c *** Path routine occasionally repeats points.
c *** Identify duplicate points and remove them!
	k=1
	raydtmp(1)=theta_path_deg(1)
	ztmp(1)=rad_path(1)
	do i=2,ntheta
	   if(abs(theta_path_deg(i)-theta_path_deg(i-1)).gt.toler)then
	      k=k+1
	      raydtmp(k)=theta_path_deg(i)
	      ztmp(k)=rad_path(i)
	   endif
	enddo
	ntheta=k
	do i=1,ntheta
	   theta_path_deg(i)=raydtmp(i)
	   rad_path(i)=ztmp(i)
	enddo

c *** find points on great circle path and stores them in xray and yray
c *** xray and yray are arrays with same length as theta_patth_deg (=ntheta)
c *** Note a slight inaccuracy as station not included in *_path arrays


c
c *** The great circle routine gceq will fail when eqlat and eqlon are EXACTLY
c *** the same. Catch this error and shift eq loc by a very tiny amount
c     

c
c        if (eqlat.eq.stlat) then
c           print*,"WARNING: eq lat was slightly shifted"
c           eqlat=eqlat-0.001
c        elseif (eqlon.eq.stlon) then
c           print*,"WARNING: eq lon was slightly shifted"
c           eqlat=eqlat-0.001
c        end if

c        print*,eqlon,eqlat
c        print*,stlon,stlat
c

        ierr=0
	call gceq(eqlon,eqlat,stlon,stlat,theta_path_deg,ntheta,xray,yray,ndim) 
        if(.not.(abs(xray(1).gt.0.d5)))then
           print*,"Datum discarded due to error in GCEQ"
           ierr=1
        end if
        if(ierr.ne.0) return
        

c        print*,xray(1),yray(1)

cTEST--plot-raypath-to-check-if-projection-worked
c                     open(771,file='./TMP/ray_test.txt')
c                     do j=1,ntheta
c                        thpath(j)=xray(j)
c                        radpath(j)=yray(j)
c                        write(771,"(f20.7,1x,5(f14.7,1x))"),thpath(j),radpath(j)
c                     enddo
c                     close(771)
c                     stop

c *** linearly interpolate xray, yray and rad_path to make sure ray path is 
c *** sufficiently densely sampled, output: yi, xi, zi; delint
	do k=1,nint ! do it nint times (default = 3)
	   call linint(yray,yi,theta_path_deg,delint,ndim,ntheta)
	   call linint_x(xray,xi,theta_path_deg,delint,ndim,ntheta)
	   call linint(rad_path,zi,theta_path_deg,delint,ndim,ntheta)
	   do i=1,ntheta*2-1
	      if(xi(i).lt.0.)xi(i)=360.+xi(i)
	      if(xi(i).gt.360.)xi(i)=xi(i)-360.
	      xray(i)=xi(i)  ! replace old xray, yray, radpath, thetapath with interpolated ones
	      yray(i)=yi(i)
	      rad_path(i)=zi(i)
	      theta_path_deg(i)=delint(i)
	   enddo
	   ntheta=2*ntheta-1 ! double the ntheta value
	enddo

	
c *** 
	do ih=1,nsqtot(nlatzones+1)
	   isave(ih)=0
	enddo


c *** Compute length of ray dl=sqrt(dr+dth*r) and put it in array xleng
        xl=0.d0
        xleng(1)=0.d0
        do i=1,ntheta-1
	   dr0=rad_path(i+1)-rad_path(i)
	   dtheta=(theta_path_deg(i+1)-theta_path_deg(i))*pi/180.d0
	   dx0=dtheta*(rad_path(i+1)+rad_path(i))*0.5d0
	   dl=sqrt(dx0*dx0+dr0*dr0)
	   xl=xl+dl
	   xleng(i+1)=xl
        enddo
        xleng(ntheta+1)=xleng(ntheta)

c *** Cubic spline interpolation of radius and theta with respect to xleng
      call drspln(1, ntheta, theta_path_deg, rad_path, qrad, fwork) ! should maybe just do linear
      call drspln(1, ntheta, xleng, rad_path, qrad_l, fwork)
      call drspln(1, ntheta, xleng, theta_path_deg, qth, fwork)
      call drspln(1, ntheta, theta_path_deg, xleng, qth_l, fwork) 



c *** Identify in which layer our ray path starts
	do ila=1,nlay
	   if((rad_path(1).lt.rbnd(ila-1)).and.(rad_path(1).ge.rbnd(ila)))then ! bug removed 13.04. la
	   iv0=ila
	   endif
      enddo


c *** Identify horizontal index of starting voxel
      ih0=isqre(yray(1),xray(1),nsqrs,nsqtot,nlatzones,n1layer,eq_incr)
      

c *** Start from pixel ih0
      isave(ih0)=1
      ind0=(iv0-1)*n1layer+ih0
      call span(ih0,ymi0,yma0,xmi0,xma0,nsqrs,nsqtot,nlatzones,eq_incr)
      x0=xray(1)
      y0=yray(1)
      z0=rad_path(1)
      d0=0.
      irec=0
	

c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c * * * * * * *  MAIN LOOP OVER ALL POINTS ON (REFINED) RAYPATH * * * * * * * * * *
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

      do i=1,ntheta !loop over all points in ray path

	if(xma0.eq.360.)xma0=0.
	if(xmi0.eq.0.)xmi0=360.

c *** Tetermine index of voxel for i-th point on ray path.
c *** The index ind is computed from iv (= vertical index) 
c *** and ih (=horizontal index).
	   do ila=1,nlay
            if((rad_path(i).lt.rbnd(ila-1)).and.(rad_path(i).ge.rbnd(ila)))then ! bug removed 13.04. la
               iv=ila
            endif
         enddo

c *** this is a hook to make sure that we neglect sensitivity in the core
	   if(rad_path(i).lt.rbnd(nlay))then
                iv=nlay+1
           endif     
              
	   ih=isqre(yray(i),xray(i),nsqrs,nsqtot,nlatzones,n1layer,eq_incr)
	   call span(ih,ymi,yma,xmi,xma,nsqrs,nsqtot,nlatzones,eq_incr)
	   ind=(iv-1)*n1layer+ih
         isave(ih)=1


c *** Check if the current point on path (in this loop ) is located inside  
c *** another voxel than the previous point, by checking if their indices
c *** are the same (ind=ind0?). If yes -> execute the following which 
c *** serves to compute dint (=degree interval) from zint and yint

	    if(ind.ne.ind0)then ! crossed over to another voxel, for endif see (*)
	        if(iv.ne.iv0)then ! vertical intersection
	            if(rad_path(i).eq.z0)then
	               print*,"Exception vert",i,ntheta,rad_path(i),z0
	               ind=ind0
	               goto 11
	            endif
	            ivint=min(iv,iv0)
	            zint=rbnd(ivint)
		    xint=x0+(xray(i)-x0)*(zint-z0)/(rad_path(i)-z0)
		    yint=y0+(yray(i)-y0)*(zint-z0)/(rad_path(i)-z0)
		    dint=d0+(theta_path_deg(i)-d0)*(zint-z0)/(rad_path(i)-z0)
	    else ! horizontal intersection
		  if(dabs(ymi-yma0).lt.toler)then
	            if(yray(i).eq.y0)then
	               print*,"Exception s to n"
	               ind=ind0
	               goto 11
	            endif
		      yint=ymi
		      xint=x0+(xray(i)-x0)*(yint-y0)/(yray(i)-y0)
		      zint=z0+(rad_path(i)-z0)*(yint-y0)/(yray(i)-y0)
		      dint=d0+(theta_path_deg(i)-d0)*(yint-y0)/(yray(i)-y0)        
		  elseif(dabs(yma-ymi0).lt.toler)then
	            if(yray(i).eq.y0)then
	               print*,"Exception n to s"
	               ind=ind0
	               goto 11
	            endif
		      yint=yma
		      xint=x0+(xray(i)-x0)*(yint-y0)/(yray(i)-y0)
		      zint=z0+(rad_path(i)-z0)*(yint-y0)/(yray(i)-y0)
		      dint=d0+(theta_path_deg(i)-d0)*(yint-y0)/(yray(i)-y0)
		  elseif(dabs(xmi-xma0).lt.toler)then
	            if(xray(i).eq.x0)then
	               print*,"Exception e to w"
	               ind=ind0
	               goto 11
	            endif
	            xint=xmi
		      yint=y0+(yray(i)-y0)*(xint-x0)/(xray(i)-x0)
		      zint=z0+(rad_path(i)-z0)*(xint-x0)/(xray(i)-x0)
		      dint=d0+(theta_path_deg(i)-d0)*(xint-x0)/(xray(i)-x0)
		  elseif(dabs(xma-xmi0).lt.toler)then
	            if(xray(i).eq.x0)then
	               print*,"Exception w to e"
	               ind=ind0
	               goto 11
	            endif
	            xint=xma
		      yint=y0+(yray(i)-y0)*(xint-x0)/(xray(i)-x0)
		      zint=z0+(rad_path(i)-z0)*(xint-x0)/(xray(i)-x0)
		      dint=d0+(theta_path_deg(i)-d0)*(xint-x0)/(xray(i)-x0)
		  else
		      print*,"Problem at point",i
		      print*,xmi0,xma0,ymi0,yma0
		      print*,xmi,xma,ymi,yma
		      print*,dabs(xmi-xma0)
		      print*,dabs(xma-xmi0)
		      print*,dabs(yma-ymi0)
		      print*,dabs(ymi-yma0)
                      ierr=1
		  endif
	    endif
   	    if(ierr.ne.0)then
	       return
	    endif



c *** Generate vector xlseg by extracting nptint (normally 10) points from the raypath within
c *** the voxel using drsple, which evaluates (extracts) values from a spline interpolated fct.
c *** Note that here use is made of the (1D) radial symmetry of my problem. Rays are traced in
c *** a spherical earth, so it is sufficient to estimate xlseg just from xleng (which is de-
c *** termined from rad_path and theta_path.

	      dtheta=dint-d0
	      dincr=dtheta/(nptint*1.d0)
	      xlseg(1)=drsple(1,ntheta,theta_path_deg,xleng,qth_l,d0)
	      do iptint=2,nptint
	         d1=d0+iptint*dincr
	         xlseg(iptint)=drsple(1,ntheta,theta_path_deg,xleng,qth_l,d1)
	      enddo




c *** Radcontribute actually computes kernel values using calcqvec and calcdrdl
c *** Note: We are stil within the loop over raypath

	      ierr=0
	      call radcontribute_drdl_vox(xleng,xlseg,nptint,qrad_l,qth,xtu,ierr)  
	      if(ierr.ne.0)return


c *** only parameterize the mantle and ignore sensitivity to core structure
c *** this is necessary to allow the use of SKS and SKKS phases, which led
c *** to wrong indices in lapos original version of voxel_int.f

	      if(iv0.le.nlay)then

                 ! Here the rows are actually updated! arowk is obtained from radcontribute
	         irec=irec+1
	         indro(irec)=ind0

                 ! find azimuth of raypath, added by ludwig
                 xlat_new=yray(i)
                 xlon_new=xray(i)

                 if (xlon_new.lt.0.0_8) xlon_new = xlon_new + 360.0_8              
                 if (i/=ntheta) then
                    xlat_old=yray(i-1)
                    xlon_old=xray(i-1)
                 end if
                 if (xlon_old.lt.0.0_8) xlon_old = xlon_old + 360.0_8               
                 call delazs(xlat_old,xlon_old,xlat_new,xlon_new,deltax,az_old,az_new)
                 if (az_new.gt.360.0_8) az_new = az_new - 360.0_8
                 if (az_new.gt.180.0_8) az_new = az_new - 180.0_8

                 if(az_new>=0.and.az_new<45) then
                    angle_cat=1
                    angle_hitcount(angle_cat,ind0)=angle_hitcount(angle_cat,ind0)+1
                    angle_hitcount(angle_cat,ind0+nvx)=angle_hitcount(angle_cat,ind0+nvx)+1
                 elseif(az_new>=45.and.az_new<90) then
                    angle_cat=2
                    angle_hitcount(angle_cat,ind0)=angle_hitcount(angle_cat,ind0)+1
                    angle_hitcount(angle_cat,ind0+nvx)=angle_hitcount(angle_cat,ind0+nvx)+1
                 elseif(az_new>=90.and.az_new<135) then
                    angle_cat=3
                    angle_hitcount(angle_cat,ind0)=angle_hitcount(angle_cat,ind0)+1
                    angle_hitcount(angle_cat,ind0+nvx)=angle_hitcount(angle_cat,ind0+nvx)+1
                 elseif(az_new>=135.and.az_new<180) then
                    angle_cat=4
                    angle_hitcount(angle_cat,ind0)=angle_hitcount(angle_cat,ind0)+1
                    angle_hitcount(angle_cat,ind0+nvx)=angle_hitcount(angle_cat,ind0+nvx)+1
                 else	
                    print*,"angle definition didn't work:"
                    print*,"ray data:", eplat, eplon, stlat, stlon !, nsqrs, nsqtot, nlatzones, numto
                    print*,xlat_old, xlon_old, xlat_new, xlon_new, az_new, az_old
                 endif

                 ! ipar 1: vph, 2:, vpv, 3: vsh, 4: vsv
                 do ipar=1,4
                    rowps(ipar,irec)=arowk(ipar,1)
                 enddo
                 
              endif  ! -> otherwise: ignore sensitivity

	      x0=xint !update x0,y0,z0 
	      y0=yint
	      z0=zint
	      d0=dint
	      iv0=iv !update other variables
	      ih0=ih !not strictly needed
	      ind0=ind 
	      xmi0=xmi
	      xma0=xma
	      ymi0=ymi
	      yma0=yma

	   endif ! Executed if intersection see (*)
c *** End of major if: above only executed when a new voxel is reached


11	   continue
	enddo !end of loop over ray path
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 



c *** Some exceptions to deal with the receiver voxel
	if(theta_path_deg(ntheta).ne.d0)then

	   dtheta=theta_path_deg(ntheta)-d0
	   dincr=dtheta/(nptint*1.d0)
	   xlseg(1)=drsple(1,ntheta,theta_path_deg,xleng,qth_l,d0)

	   do iptint=2,nptint
	      d1=d0+iptint*dincr
	      if(d1.ge.theta_path_deg(ntheta))d1=d1-0.001 ! Not sure why there is this stupid problem <-lb
	      xlseg(iptint)=drsple(1,ntheta,theta_path_deg,xleng,qth_l,d1)
	   enddo
           ierr=0

	   call radcontribute_drdl_vox(xleng,xlseg,nptint,qrad_l,qth,xtu,ierr)
	   if(ierr.ne.0)then
	      return
	   endif

	   irec=irec+1
	   indro(irec)=ind0
           
           ! still entirely general 
           do ipar=1,4
              rowps(ipar,irec)=arowk(ipar,1)
           enddo

	endif
        return
	end


c
c check if value is nan
c
      logical function inan(x)
        real x
        inan = .true.
        If(x .ge. 0.0) then
           inan = .false.
        elseif (x .le. 0.0) then
           inan = .false.
        endif
      end function inan


c
c subroutine delazs
c
      subroutine delazs(eplat,eplong,stlat,stlong,delta,azep,azst)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      data hpi,twopi,rad,reprad/1.5707963268
     1   ,6.2831853,.017453293,57.2957795/
      arcos(x)=atan2(sqrt(1.-x*x),x)
      el=eplat*rad
      el=hpi-el
      stl=stlat*rad
      stl=hpi-stl
      elon=eplong*rad
      slon=stlong*rad
      as=cos(stl)
      bs=sin(stl)
      cs=cos(slon)
      ds=sin(slon)
      a=cos(el)
      b=sin(el)
      c=cos(elon)
      d=sin(elon)
      co=1.0
      cdel=a*as+b*bs*(c*cs+d*ds)
      if(abs(cdel).gt.1.0) cdel=sign(co,cdel)
      delt=arcos(cdel)
      delta=delt*reprad !distance
      sdel=sin(delt)
      caze=(as-a*cdel)/(sdel*b)
      if(abs(caze).gt.1.0) caze=sign(co,caze)
      aze=arcos(caze)
      if(bs.gt.0.0) cazs=(a-as*cdel)/(bs*sdel)
      if(bs.eq.0.0) cazs=sign(co,cazs)
      if(abs(cazs).gt.1.0) cazs=sign(co,cazs)
      azs=arcos(cazs)
      dif=ds*c-cs*d
      if(dif.lt.0.0) aze=twopi-aze
      azep=reprad*aze
      if(dif.gt.0.0) azs=twopi-azs
      azst=reprad*azs
      return
      end
