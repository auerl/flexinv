      subroutine addsurfspline(eqlat,eqlon,stlat,stlon,
     +	bslat,bslon,ngpt,h,arytrans,ntrans,aryrefl,nrefl,
     +	pparm,rname,iray,rdisp,indskip)
      implicit double precision (a-h,o-z)
c  This subroutine locates and computes surface kernels for 
c  a reflected phase such as SdS or PdP at a discontinuity.
c
c	 eqlat ---  Earthquake source lat
c	 eqlon ---  Earthquake source lon
c	 stlat ---  Earthquake receiver lat
c	 stlon ---  Earthquake receiver lon
c	 bslat ---  lat array of horizontal splines
c	 bslon ---  lon array of horizontal splines
c	 ngpt  ---  number of horizontal grid points
c	 h     ---  bucket size (arc length of nonzero B-spline contribution)
      parameter(maxfreq=8)
      parameter(maxnode=maxfreq*maxfreq*10+2)
      parameter(maxrknot=50)
      parameter(maxparm=maxnode*maxrknot)
      data pi /3.141592653579/
      data rad/57.2957795130824/
      data reprad /57.295779579/    


      common/partvel/pathder(5,3200),pathvel(5,3200)
      common/amatrix/arow(maxparm),indarow(maxparm),nonzero_a
      common/invopt/numinvvar,invparlist(5)
      common/radcontrib/qseg(5,3200),arowk(5,maxrknot)
      common/anisopath/theta_path(3200),rad_path(3200),ntheta
      common/layr/nl(20),xb(20),xt(20),ifanis,nplay
c***  secondary phases ------------------------------------------------
c
      common/secphases/isectype,rsec,indexsec,lsec
c----------------------------------------------------------------------
      dimension bslat(1),bslon(1),arytrans(1),aryrefl(1)
      dimension transmat(9)
      character*(*) rname

      hpi=pi*0.5d0
      twoh=2.d0*h
      twohsq=twoh*twoh
      call conv2geocen(eqlat, eqlon, eqtheta, eqphi)
      call conv2geocen(stlat, stlon, sttheta, stphi)
      call calc_fwdtransmat(eqtheta, eqphi, sttheta, stphi, transmat)
      bstrans=0.d0
      bstrans1=0.d0
      bsrefl=0.d0
      if(ntrans.ne.0) bstrans=surftransker(pparm,rdisp,rname,bstrans1,ibmode)
      if(nrefl.ne.0) bsrefl=surfreflker(pparm,rdisp,rname)
c	if(ntrans.ne.0) print*, 'bstrans=', bstrans
c        if(nrefl.ne.0) print*, '  bsrefl=', bsrefl
c	write(113, *) eqlon, eqlat
c	write(113, *) stlon, stlat
      nbs_r=0
      nbs_t=0
      do i=1, ngpt
	ii=indskip+i
        call conv2geocen(bslat(i), bslon(i), sptheta, spphi)
c... transform spline location to source--receiver coordinates
        call transform(sptheta, spphi, transmat, spth1, spph1)
	if(abs(spth1-hpi).le.twoh) then
	   aa=hpi-spth1
	   do j=1, nrefl
c... account for reflection at discontinuity
	   	bb=aryrefl(j)-spph1
	   	dist=sqrt(aa*aa+bb*bb) ! distance of bell to reflection pt
		if(dist.le.twoh) then
			nbs_r=nbs_r+1
			bsp=fdelta(dist,h)
			if(iray.eq.1) then
				arow(ii)=arow(ii)+bsp*bsrefl
			else
				arow(ii)=arow(ii)-bsp*bsrefl
			endif
c			write(114, *) bslon(i), bslat(i), bsp
		endif
	   enddo
	   do j=1, ntrans
c... account for transmission at discontinuity
	   	bb=arytrans(j)-spph1
	   	dist1=sqrt(aa*aa+bb*bb) ! distance of bell to reflection pt
		if(dist1.le.twoh) then
			nbs_t=nbs_t+1
			bsp1=fdelta(dist1,h)
			if(ibmode.eq.0) then
c... regular phase such as S, ScS, or for SS or PP precursors
			   xker=bsp1*bstrans
			else if(ibmode.eq.1) then
c... PS or SP, second leg is different from first
			   if(j.le.2) then
				xker=bsp1*bstrans
			   else 
				xker=bsp1*bstrans1
			   endif
			else
c... PdS or SdP, second point is the conversion point
			   if(j.eq.1) then
				xker=bsp1*bstrans
			   else 
				xker=bsp1*bstrans1
			   endif
			endif
			if(iray.eq.1) then
				arow(ii)=arow(ii)+xker
			else
				arow(ii)=arow(ii)-xker
			endif
c			write(115, *) bslon(i), bslat(i), bsp1
		endif
	   enddo
	endif
      enddo
      return
      end
