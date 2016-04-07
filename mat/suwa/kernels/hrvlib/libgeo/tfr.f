C     *** TFR ********* VERSION 1, MODIFICATION LEVEL 0 *** DKO20813 ***
C     *                                                                *
C     *   COMPLETE REAL FOURIER ANALYSIS IF IOPT IS NEGATIVE           *
C     *   INITIALIZE REAL FOURIER SYNTHESIS IF IOPT IS POSITIVE        *
C     *                                                                *
C     *   5736-XM7 COPYRIGHT IBM CORP. 1971                            *
C     *   REFER TO INSTRUCTIONS ON COPYRIGHT NOTICE FORM NO. 120-2083  *
C     *   FE SERVICE NO. 200281                                        *
C     *                                                                *
C     ******************************************************************
C
      SUBROUTINE TFR(RA,RB,N,IOPT,IER)
      DIMENSION RA(1),RB(1)
      DOUBLE PRECISION CD,DA,DI,DN,DV
c      EQUIVALENCE (DI,CI),(DN,SN)
C
      ISW  =IER
      IER  =0
      IF (N-1) 1,1,2
    1 IER  =1000
      GOTO 4
    2 IF (IOPT) 6,3,6
    3 IER  =2000
    4 IF (ISW+12345) 5,14,5
    5 CALL WIER(IER,20813)
      GOTO 14
    6 NJ   =IABS(IOPT)
      L    =N*NJ+1
      DN   =1.D0/DFLOAT(N)
      sn=dn
      DA   =3.1415926535897932D0*DN
      AA   =RA(1)
      IF (IOPT) 7,7,10
    7 BB   =RB(1)
      RA(1)=SN*(AA+BB)
      RA(L)=SN*(AA-BB)
      K    =N/2
      IF (N-K-K) 8,8,9
    8 K    =K*NJ+1
      RA(K)=SN*RA(K)
      RB(K)=SN*RB(K)
    9 RB(1)=0.
      RB(L)=0.
      DN   =0.5D0*DN
      sn=dn
      GOTO 11
   10 BB   =RA(L)
      DN   =0.5D0
      sn=dn
      RA(1)=SN*(AA+BB)
      RB(1)=SN*(AA-BB)
   11 IF (N-2) 14,14,12
   12 CD   =DCOS(DA)
      SD   =DSIN(DA)
      SC   =DN
      DI   =DN*CD
      ci=di
      SI   =SN*SD
      DV   =DI
      DA   =DN
      K    =L
      J    =NJ+1
      NH   =((N-1)/2)*NJ+1
      DO 13 I=J,NH,NJ
         K    =K-NJ
         AA   =RA(I)
         AB   =RA(K)
         AC   =SC*(AA+AB)
         AB   =AA-AB
         BA   =RB(I)
         BB   =RB(K)
         BC   =BA+BB
         BA   =SC*(BA-BB)
         AA   =CI*BC-SI*AB
         BB   =CI*AB+SI*BC
         RA(I)=AC+AA
         RA(K)=AC-AA
         RB(I)=BB-BA
         RB(K)=BB+BA
         DN   =DV+DI
         sn=dn
         SI   =SN*SD
         DV   =DN*CD
         DI   =DV-DA
         ci=di
   13    DA   =DN
   14 RETURN
      END
