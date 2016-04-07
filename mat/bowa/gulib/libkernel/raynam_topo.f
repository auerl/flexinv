      subroutine raynam_topo(ray,ierr)
      implicit real*8(a-h,o-z)
c
c only names using P,S,I,J,K,p,s,i,c allowed
c p, P, s, or S must be first character.
c The reflection and turning ray , e.g. P and PcP are distinguished 
c by the ray parameter range.
c This version handles the SS precursors and converted waves off of
c 220, 400, 670 discontinuities, JG2002.
c
c
      common/ray$/pmin,pmax,nps(20,2),npsrc(2),nn,iup,jref,npleg,nsleg
      common/slow/picblp,picbls,picbup,pcmblp,pcmbup,pcmbus,
     +   p660lp,p660ls,lic,loc,llm,lsrc,klay(20)
      common/layr/nl(20),xb(20),xt(20),ifanis,nplay
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
c
c secondary phases
c
      common/secphases/isectype,rsec,indexsec,lsec

      logical split
      character*(*) ray      
      data ptol/1.d-6/

c
c split is false if source is at existing interface (eg surface)
c if split is true, we put multiplicity for legs above source in
c npsrc, rdep = source radius
c xt is the terminating radius for the shells in polyprem
c

      split=.true.
      if(rdep.eq.xt(lsrc)) split=.false.

c
c get indices for bottom of layers -- note
c that the source layer and the turning layer will 
c have two and one extra points respectively.
c

      klay(1)=nl(1)-1
c
c the following is problematic part of the old code, Yu Gu
c      if(lsrc.eq.1.and.split) klay(1)=klay(1)+2
c      if(lsrc.eq.1.and.split) then
c		print*, 'nl(1)=', nl(1), ' klay=', klay(1)+2, '   klay corrected=', klay(1)*2
c		klay(1)=klay(1)*2
c      endif
c
      do i=2,numlyr
         klay(i)=nl(i)+klay(i-1)-1
c
c the following is problematic part of the old code, Yu Gu
c        if(lsrc.eq.i.and.split) klay(i)=klay(i)+2
c
         if(lsrc.eq.i.and.split) then
		klay(i)=(nl(i)*2)+klay(i-1)-2
	 endif
      enddo
      ierr=0                                
c
c iup=1 if little p or 2 if little s
c
      iup=0

c
c jref=1 means reflection from cmb, 2= reflection from icb
c
      jref=0
      npsrc(1)=0
      npsrc(2)=0
      do i=1,numlyr
        nps(i,1)=0
        nps(i,2)=0  
      enddo                                                  


c ************************************************* evaluate # ray legs
      j1=0
      do 10 i=1,70
c
c after ichdec is called, if j=11, then no matching phase found
c if j=1, then empty name is found.
c otherwise, ' PSIJKpsic' determines the return value j.
c k = ray type (1=compression wave,  2=shear wave), JG2002.
c

      j=ichdec(ray(i:i),k)
      goto(10,11,11,12,12,13,14,14,15,16,20),j
c**** P or S ****
   11 goto(20,111,112,20,20,114,113,113,20,100),j1
c**** P or S (only done for beginning segment) ****
      call addray(nps(1,k),lsrc,loc-1,1)
      goto 100
c------------------------------------------------------------------------
c**** previous segment is P wave
  111 if(isectype.eq.0) then
c**** if a major phase such as PS or PP
      	call addray(nps(1,1),1,loc-1,1)
      	if(split) npsrc(1)=npsrc(1)+1
      	call addray(nps(1,k),1,loc-1,1)
      	if(split) npsrc(k)=npsrc(k)+1
      else
      	call addray(nps(1,1),lsec,loc-1,1)
      	if(split.and.rdep.le.rsec) npsrc(1)=npsrc(1)+1
	if(isectype.eq.1) then
c**** PP precursors
      		call addray(nps(1,k),lsec,loc-1,1)
      		if(split.and.rdep.le.rsec) npsrc(k)=npsrc(k)+1
      	else if(isectype.eq.2) then
c**** if PS conversions, S segment
      		call addray(nps(1,k),1,lsec-1,1)
      		if(split.and.rdep.gt.rsec) npsrc(k)=npsrc(k)+1
	endif
      endif         
      goto 100
c------------------------------------------------------------------------
c**** previous segment is S wave
  112 if(isectype.eq.0) then
c**** major phase
	call addray(nps(1,2),1,loc-1,1)         
      	if(split) npsrc(2)=npsrc(2)+1
      else
c**** SS precursors or SP conversions
	call addray(nps(1,2),lsec,loc-1,1)         
      	if(split.and.rdep.le.rsec) npsrc(2)=npsrc(2)+1
      endif
  113 if(isectype.eq.0) then
c**** major phase
	call addray(nps(1,k),1,loc-1,1)              
      	if(split) npsrc(k)=npsrc(k)+1
      else
	if(isectype.eq.1) then
c**** SS precursors
		call addray(nps(1,k),lsec,loc-1,1)              
      		if(split.and.rdep.le.rsec) npsrc(k)=npsrc(k)+1
	else if(isectype.eq.2) then
c**** SP conversions, P segment
		call addray(nps(1,k),1,lsec-1,1)              
      		if(split.and.rdep.gt.rsec) npsrc(k)=npsrc(k)+1
	endif
      endif
      goto 100

c------------------------------------------------------------------------
  114 call addray(nps(1,1),loc,lic-1,1)
      goto 100
c------------------------------------------------------------------------

C**** I or J ****
   12 goto(20,20,20,121,121,121,20,20,20,20),j1
  121 call addray(nps(1,k),lic,numlyr,2)           
      goto 100
C**** K ****
   13 goto(20,132,132,100,100,131,20,20,100,20),j1
  131 call addray(nps(1,1),loc,lic-1,1)            
  132 call addray(nps(1,1),loc,lic-1,1)            
      goto 100
C**** p or s ****
   14 if(j1.ne.0) goto 20
      iup=k
      call addray(nps(1,k),1,lsrc-1,1)
      if(split) npsrc(k)=npsrc(k)+1
      goto 100
c**** i ****
   15 jref=2
      if(j1.ne.6) goto 20
      goto 100
C**** c ****
   16 jref=1
      if(j1.lt.2.or.j1.gt.3) goto 20
  100 j1=j  
      k1=k
   10 continue
c----------------------------------------------------------------
C**** receiver leg ****
      goto(20,201,201,20,20,20,201,201,20,20),j1
  201 if(isectype.ne.2) then
c**** if not converted phases, add final segment
      	call addray (nps(1,k1),1,loc-1,1)
      	if(split) npsrc(k1)=npsrc(k1)+1
      endif
c----------------------------------------------------------------
c*****finish of ray leg summation
c**** nn is total # layers sampled
      do j=numlyr,1,-1
        nn=j
        if(nps(j,1)+nps(j,2).ne.0) goto 300
      enddo
  300 npleg=0
      nsleg=0
      do i=1,nn
        npleg=npleg+nps(i,1)
        nsleg=nsleg+nps(i,2)
      enddo
c**** evaluate ray parameter range
      if(jref.eq.1.and.nn.lt.loc) then
c**** bounce off core (no legs beneath unlike e.g. PcPPKP)
        pmax1=1.e10
        pmax2=1.e10
        if(nps(loc-1,1).ne.0) pmax1=pcmbup
        if(nps(loc-1,2).ne.0) pmax2=pcmbus
        pmax=min(pmax1,pmax2)-ptol
        pmin=0.
c**** bounce off inner core -- no legs beneath
      else if(jref.eq.2.and.nn.lt.lic) then
        pmax=picbup-ptol
        pmin=0.
      else
c**** turning rays
c**** following for rays turning in outer core
        if(nn.eq.loc) then
          pmin=picbup+ptol
          pmax1=1.e10
          pmax2=1.e10
          if(nps(loc-1,2).ne.0) pmax1=pcmbus
          if(nps(loc-1,1).ne.0) pmax1=pcmbup
          pmax=min(pmax1,pmax2)
          pmax=min(pmax,pcmblp)-ptol
          return
c**** following for rays turning in inner core
c**** program has numerical problems if we turn very close
c**** to center of Earth ... hence limit on pmin
        else if(nn.eq.lic) then
          pmin=.0001
          pmax1=1.e10
          pmax2=1.e10
          if(nps(lic,1).ne.0) pmax1=picblp
          if(nps(lic,2).ne.0) pmax2=picbls
          pmax=min(pmax1,pmax2)
          pmax=min(pmax,picbup)-ptol
          return
        else
c**** this for mantle-turning rays
          pmin1=0
          pmin2=0
          if(nps(loc-1,2).ne.0) pmin1=pcmbus
          if(nps(loc-1,1).ne.0) pmin2=pcmbup
          pmin=max(pmin1,pmin2)+ptol

c          if(iup.ne.0)then
c  p or s at source
c            y=rdep/rnorm
c            call getmod(y,lsrc,rho,vpv,vsv,vph,vsh,eta)
c            if(iup.eq.1) pmax=rdep/min(vpv,vph)
c            if(iup.eq.2) pmax=rdep/min(vsv,vsh)
c          else 
c p or s at source
c make lower-mantle turning only
c
            pmax1=1.e10
            pmax2=1.e10
            if(nps(lsrc,1).ne.0) pmax1=p660lp
            if(nps(lsrc,2).ne.0) pmax2=p660ls
            pmax=min(pmax1,pmax2)-ptol
c          endif
        end if
      end if
      return                      
C**** error in ray name
   20 ierr=1
      return
      END
