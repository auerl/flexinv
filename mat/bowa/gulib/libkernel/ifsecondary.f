      function ifsecondary(phasein, phaseout)
      implicit real*8(a-h,o-z)
c
c  Checks if the ray name is a secondary phase including 
c  SdS (SS precursors), SdP or PdS (converted phases), JG 2002.
c
      common/coeff/coef(4,8,20)
      common/layr/nl(20),xb(20),xt(20),ifanis,nplay
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      common/slow/picblp,picbls,picbup,pcmblp,pcmbup,pcmbus,
     +   p660lp,p660ls,lic,loc,llm,lsrc,klay(20)
      data pi/3.14159265358979d0/,rtol/0.d0/
c
c***  secondary phases
c
      common/secphases/isectype,rsec,indexsec,lsec
      common/slow1/l220,l400,l520
c*** isectype:  0=major phase,  1=reflected, 2=converted
c*** depsec:	depth of interface for conversion or reflection
c*** indexsec:	index of the secondary phase
c*** lnsec:	layer number for the bottom of conversion
c
      parameter(maxphase=30)
      parameter(nphase=16)
      character*(*) phasein, phaseout
      character*20 phasename(maxphase)
      data (phasename(i), i=1, nphase)/
     +  'S220S'
     + ,'S400S'
     + ,'S520S'
     + ,'S670S'
     + ,'P220P'
     + ,'P400P'
     + ,'P520P'
     + ,'P670P'
     + ,'S220P'
     + ,'S400P'
     + ,'S520P'
     + ,'S670P'
     + ,'P220S'
     + ,'P400S'
     + ,'P520S'
     + ,'P670S' /
c
c
      ifsecondary=0
      isectype=0
      phaseout=phasein
      do i=1, nphase
	if(phasein.eq.phasename(i)) then
        	ifsecondary=i
c		secray=phasein ! commented out by lapo 22.1.2010
		goto 10
	endif
      enddo
c*** if not secondary phase
c      write(*, "('Input phase = major phase')")
      return
c
c*** if secondary phase
10    indexsec=i
      if(i.lt.9) then
	isectype=1
      else 
	isectype=2
      endif

c      write(*, "('Input phase = secondary phase')")
c      write(*, "('phase number= ', i2)") indexsec
c      print*, '  name= ', phasename(i)
      do k=1,numlyr
        if(abs(xt(k)-5820.).lt.20.) l520=k
        if(abs(xt(k)-5970.).lt.20.) l400=k
        if(abs(xt(k)-6150.).lt.20.) l220=k	
      enddo
      if(i.eq.1.or.i.eq.5.or.i.eq.9.or.i.eq.13)  then
	rsec=6151.
	lsec=l220
      endif
      if(i.eq.2.or.i.eq.6.or.i.eq.10.or.i.eq.14) then
	rsec=5971.
	lsec=l400
      endif
      if(i.eq.3.or.i.eq.7.or.i.eq.11.or.i.eq.15) then
	rsec=5821.
	lsec=l400
      endif
      if(i.eq.4.or.i.eq.8.or.i.eq.12.or.i.eq.16) then
	rsec=5701.0
	lsec=llm
      endif
      if(i.le.4) then
	phaseout='SS'
      else if(i.gt.4.and.i.le.8) then
	phaseout='PP'
      else if(i.gt.8.and.i.le.12) then
      	phaseout='SP'
      else
      	phaseout='PS'
      endif
      return
      end
