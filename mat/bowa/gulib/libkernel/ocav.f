CPROG OCAV
CXREF

	SUBROUTINE ocav(nsqrs,nsqtot,nlatzones,eq_incr,n,
     &	eplat,eplon,stlat,stlon,ientm,jsqrm,ddelm,
     &	ient,jsqrg,ddelg,pi,DEL,xlatm,xlonm,xlatg,xlong)
      parameter (ndim=10000)
	dimension geo_colat(ndim),geo_long(ndim)
	DIMENSION PHT(ndim),IDX(ndim),jsqrg(ndim),ddelg(ndim)
	DIMENSION jsqrm(ndim),ddelm(ndim)
	dimension nsqrs(nlatzones),nsqtot(nlatzones)
        real*8,dimension(ndim) :: xlonm,xlong, xlatm,xlatg
        real*8 eplat, eplon, stlat, stlon     


	PI2=PI*2.
	FORPI=4.*PI
	RADIAN=180./PI
	nth1=nlatzones+1
	dth=pi/(nlatzones*1.)

      TH1=(90.-EPLAT)/RADIAN
      PH1=EPLON/RADIAN
      TH2=(90.-STLAT)/RADIAN
      PH2=STLON/RADIAN
      STH1=SIN(TH1)
      STH2=SIN(TH2)
      CTH1=COS(TH1)
      CTH2=COS(TH2)
      SPH1=SIN(PH1)
      SPH2=SIN(PH2)
      CPH1=COS(PH1)
      CPH2=COS(PH2)
      CPH21= CPH1*CPH2+SPH1*SPH2
      SPH21=SPH2*CPH1-SPH1*CPH2
      CDEL=STH1*STH2*CPH21+CTH1*CTH2
      CCAPTH=STH1*STH2*SPH21/SQRT(1.-CDEL*CDEL)
      SCAPTH=SQRT(1.-CCAPTH*CCAPTH)
      capth=atan2(scapth,ccapth)
      SCAPPH=CTH1*STH2*CPH2-CTH2*STH1*CPH1
      CCAPPH=STH1*CTH2*SPH1-STH2*CTH1*SPH2
      CAPPH=ATAN2(SCAPPH,CCAPPH)
      SCAPPH=SIN(CAPPH)
      CCAPPH=COS(CAPPH)
C
      DEL=ATAN2(SQRT(1.-CDEL*CDEL),CDEL)
      CPHSP=CCAPTH*STH1*(CPH1*CCAPPH+SPH1*SCAPPH)-SCAPTH*CTH1
      SPHSP=STH1*(SPH1*CCAPPH-CPH1*SCAPPH)
      PHSP=ATAN2(SPHSP,CPHSP)
c      write(*,*) del*radian, capth*radian, capph*radian,
c     1phsp*radian
      thet=capth
      if(capth.gt.0.5*pi) thet=pi-capth
      thmin=0.5*pi-thet
      thmax=0.5*pi+thet
      lat_zone=0
      IENT=0
      IF(SCAPTH.EQ.0) GOTO 10
C***
      DO 20 I=2,NTH1
      lat_zone=lat_zone+1
      TH=FLOAT(I-1)*DTH
      CTH=COS(TH)
      CPHT=-CTH/SCAPTH
      CPHT2=CPHT*CPHT
      if(th.lt.thmin.or.th-dth.gt.thmax) then
         go to 20
      endif
      IF(CPHT2.GT.1.) GOTO 21
      IENT=1+IENT
      PHT(IENT)=ATAN2(SQRT(1.-CPHT2),CPHT)
      spht=sin(pht(ient))
      geo_colat(ient)=th
      geo_long(ient)=atan2(spht,-cth*ccapth/scapth)
      xel=AMOD(PHT(ient)-phsp+FORPI,PI2)
      IENT=1+IENT
      PHT(IENT)=-PHT(IENT-1)
      geo_colat(ient)=th
      geo_long(ient)=atan2(-spht,-cth*ccapth/scapth)
      if(geo_long(ient).lt.0) geo_long(ient)=pi2+geo_long(ient)
      xel=AMOD(PHT(ient)-phsp+FORPI,PI2)
  21  numlong=nsqrs(lat_zone)
      dphi=pi2/float(numlong)
      DO 40 j=1,numlong
c	print*,j,numlong
c	if((j.gt.numlong).or.(j.lt.1))pause
         PH=FLOAT(j-1)*DPHi
         angr=ph-capph
         thlo=atan(-ccapth/(scapth*cos(angr)))
         if(thlo.lt.0.) thlo=pi+thlo
         if(thlo.gt.th-dth.and.thlo.lt.thmin) go to 40
         if(thlo.lt.th+dth.and.thlo.gt.thmax) go to 40
         if(thlo.gt.th.or.thlo.lt.th-dth) go to 40
         SPH=SIN(PH)
         CPH=COS(PH)
         IENT=IENT+1
         PHT(IENT)=ATAN2(CCAPTH*(SPH*CCAPPH-CPH*SCAPPH),
     1        CPH*CCAPPH+SPH*SCAPPH)
         IF(PHT(IENT).GT.PI) PHT(IENT)=PHT(IENT)-PI2
         geo_colat(ient)=thlo
         geo_long(ient)=ph
         xel=AMOD(PHT(ient)-phsp+FORPI,PI2)
   40 CONTINUE
   20 CONTINUE

   10 CONTINUE
C

!       PRINT*, IENT
!       PAUSE
 
      DO 60 I=1,IENT
   60 PHT(I)=AMOD(PHT(I)-PHSP+FORPI,PI2)
      IENT=IENT+1
      PHT(IENT)=0.
      IENT=IENT+1
      PHT(IENT)=DEL
C
      CALL RSOINC(PHT,IENT,IDX)
C

c      write(88,*) pht
      

      PHT(IENT+1)=PI2
      ientm=0
      ientg=0
      cdelta=0.
      DO 50 I=1,IENT
      I1=1+I
      PHTT=.5*(PHT(I)+PHT(I1))+PHSP
      CPHT=COS(PHTT)
      SPHT=SIN(PHTT)
      CTH=-CPHT*SCAPTH
      TH=ATAN2(SQRT(1.-CTH*CTH),CTH)
      CPH=CPHT*CCAPTH*CCAPPH-SPHT*SCAPPH
      SPH=CPHT*CCAPTH*SCAPPH+SPHT*CCAPPH
      IF(SPH.EQ.0..AND.CPH.EQ.0.) GOTO 9873
      PH=ATAN2(SPH,CPH)
      if(ph.lt.0.) ph=ph+pi2
      if(ph.gt.pi2) ph=ph-pi2
      GOTO 9874
 9873 PH=0.
 9874 CONTINUE
      XLAT=90.-TH*RADIAN
      ilat=xlat*100.+0.5
      XLON=PH*RADIAN
      ilon=xlon*100.+0.5
c      jsqre=SUPERisqre(ilat,ilon,nsqrs,nsqtot,nlatzones,n,eq_incr)
	jsqre=isqre(xlat,xlon,nsqrs,nsqtot,nlatzones,n,eq_incr) ! lapo 16/5/2006
      RD=(PHT(I1)-PHT(I))
      IF(PHT(I1).LE.DEL) then
      ientm=ientm+1
      jsqrm(ientm)=jsqre
      ddelm(ientm)=rd
      xlatm(ientm)=xlat
      xlonm(ientm)=xlon
c--correction
      if(ientm.gt.1.and.jsqrm(ientm).eq.jsqrm(ientm-1)) then
c	print*,'repeated index',jsqrm(ientm)
      ientm=ientm-1
      ddelm(ientm)=ddelm(ientm)+rd
      endif
c------------
      cdelta=cdelta+rd
	call range(jsqre,XLAMIN,XLAMAX,blocla,
     &	XLOMIN,XLOMAX,bloclo,nsqrs,nsqtot,nlatzones,n,eq_incr)
      endif

      ientg=ientg+1
      jsqrg(ientg)=jsqre
      ddelg(ientg)=rd
      xlong(ientg)=xlon
      xlatg(ientg)=xlat
      if(ientg.gt.1.and.jsqrg(ientg).eq.jsqrg(ientg-1)) then
      ientg=ientg-1
      ddelg(ientg)=rd+ddelg(ientg)
      endif
   50 continue
      RETURN
      END
