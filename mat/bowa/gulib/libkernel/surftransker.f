      double precision function surftransker(rparm,rdis,rname,surfaux,ibmode)
      implicit double precision (a-h,o-z)
c..  dt due to dr Morrelli and Dziewonski
c... note I took the negative of the above due to depth, not radius
c... surfaux takes consideration for converted phases, both for regular
c... phases like SP or PS, but for SdP and PdS.  It puts the 
c... value of one value in surftransker, but auxilary one in surfaux.  
c... They are different since PdS or SdP involves two kinds wave types
c... in the second part (upgoing) of the leg.  Note this only considers
c... isotropic (actually, equivalent) velocities.
c... ibmode = 0 (regular phase or SdS or PdP, no conversions)
c... ibmode = 1 (conversion of type SP or PS)
c... ibmode = 2 (conversion of type SdP or PdS)
      common/layr/nl(20),xb(20),xt(20),ifanis,nplay
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      common/coeff/coef(4,8,20)
      common/premvel/vphp,vpvp,vshp,vsvp,etap,rhop,vpiso,vsiso
c***  secondary phases ------------------------------------------------
c
      common/secphases/isectype,rsec,indexsec,lsec
c----------------------------------------------------------------------
c
      character*(*) rname
c
      r=rdis-0.1
      y=r/rnorm
      do j=1, numlyr
	      if(r.gt.xb(j).and.r.le.xt(j)) then
		      iq = j
		      ierror = 0
		      goto 999
	      endif
      enddo
999   call getmod_iso(y,iq,vp,vs)
      r=rdis+0.1
      y=r/rnorm
      do j=1, numlyr
	      if(r.gt.xb(j).and.r.le.xt(j)) then
		      iq = j
		      ierror = 0
		      goto 1000
	      endif
      enddo
1000  call getmod_iso(y,iq,vp1,vs1)
c	print*, 'vp=', vp, '  vp1=', vp1
c	print*, 'vs=', vs, '  vs1=', vs1
c------------------------------------------------------------------------
      j=ichdec(rname(1:1),k)
      ibmode=0
      if(npleg.ne.0.and.nsleg.ne.0) then
c... converted phase
	print*, 'converted phase, RIGHT!!!'
      	if(isectype.ne.2) then
c... PS or SP
	   if((npleg+nsleg).gt.2) stop 'Cannot yet deal with phases similar to PSP, SPS, etc..!'
	   ibmode=1
	   if(k.eq.1) then
c..  PS
		vtop=vp1
		vbot=vp
		vtop1=vs1
		vbot1=vs
	   else
c..  SP
		vtop=vs1
		vbot=vs
		vtop1=vp1
		vbot1=vp
	   endif
      	else
c... PdS or SdP
	   ibmode=2
	   if(k.eq.1) then
c... PdS
		vtop=vp1
		vbot=vp
		vtop1=vs1
		vbot1=vp
	   else
c... SdP
		vtop=vs1
		vbot=vs
		vtop1=vp1
		vbot1=vs
	   endif
      	endif
      else
	if(k.eq.1) then
		vtop=vp1
		vbot=vp
	else
		vtop=vs1
		vbot=vs
	endif
      endif
c------------------------------------------------------------------------
      etatop=rdis/vtop
      etabot=rdis/vbot
      SurfTransKer=(sqrt(etatop*etatop-rparm*rparm)
     &             -sqrt(etabot*etabot-rparm*rparm))/rdis
      if(npleg.ne.0.and.nsleg.ne.0) then
c... if converted phase, need two kernels SurfTransKer and surfaux
      	 etatop=rdis/vtop1
      	 etabot=rdis/vbot1
	 print*, 'Converted:', ' vtop1 =', vtop1, '  vbot1=', vbot1
         surfaux=(sqrt(etatop*etatop-rparm*rparm)
     &             -sqrt(etabot*etabot-rparm*rparm))/rdis
      endif
      return
      end
