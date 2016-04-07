
CPROG BHOUSE
      SUBROUTINE BHOUSE(N,C,Y,A,B,P,TA,TB,W,V,EV,UU,IDU)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C   SCM VERSION
      DIMENSION A(1),B(1),P(1),TA(1),TB(1),W(1),Y(1),V(1)
     1   ,UU(IDU,IDU)
      DIMENSION C(3),EV(1)
      EPS=1.D-14
      UMEPS=1.D0-EPS
      TOL=1.D-70
      JSKIP=0
      KSKIP=1
      NM1=N-1
      I=1
      IDM1=0
      P(1)=0.D0
      V(1)=0.D0
      W(1)=0.D0
      IF(N.LE.0) RETURN
      IF(N.GT.2) GO TO 4
      IF(N.EQ.2) GO TO 3
      EV(1)=C(1)
      Y(1)=1.D0
C     WRITE(IU) EV(1),Y(1)
C     END FILE IU
C     REWIND IU
      RETURN
    3 A(1)=C(1)
      B(1)=-C(2)
      CKJ=C(3)
      IP1=2
      GO TO 215
    4 IP1=I+1
      NMI=N-I
      KJ=IDM1
      J=I
    5 JP1=J+1
      VJ=V(J)
      K=J
      LJ=N-J+1
      JD=KJ+1
      IF(KSKIP.EQ.1) GO TO 6
      PJ=P(J)
      WJ=W(J)
    6 KJ=KJ+1
      CKJ=C(KJ)
      IF(KSKIP.EQ.1) GO TO 7
      DC=-(PJ*W(K)+WJ*P(K))
      CKJ=DC+CKJ
      C(KJ)=CKJ
    7 IF(J.GT.I) GO TO 14
      IF(K.GT.J) GO TO 8
      A(I)=CKJ
      K=K+1
      GO TO 6
    8 Y(K)=0.D0
      V(K)=CKJ
      K=K+1
      IF(K.LE.N) GO TO 6
      JSKIP=0
      SUM=DOT(V(JP1),1,V(JP1),1,LJ-1)
      IF(SUM.LE.TOL) GO TO 10
      S=DSQRT(SUM)
      CSD=V(JP1)
      IF(CSD.LT.0.D0) S=-S
      V(JP1)=CSD+S
      C(JD+1)=V(JP1)
      H=SUM+CSD*S
      B(I)=-S
      GO TO 12
   10 B(I)=0.D0
      JSKIP=1
   12 IDM1=KJ
      IF(JSKIP.EQ.1.AND.KSKIP.EQ.1) GO TO 215
      J=JP1
      GO TO 5
   14 IF(JSKIP.EQ.0) GO TO 15
      K=K+1
      IF(K.LE.N) GO TO 6
      J=JP1
      IF(J.LE.N) GO TO 5
      GO TO 215
   15 Y(K)=Y(K)+CKJ*VJ
      K=K+1
      IF(K.LE.N) GO TO 6
      IF(J.EQ.N) GO TO 17
      Y(J)=Y(J)+DOT(C(JD+1),1,V(JP1),1,LJ-1)
      J=JP1
      GO TO 5
   17 SP=DOT(V(IP1),1,Y(IP1),1,NMI)/(H+H)
      DO 21 J=IP1,N
      W(J)=V(J)
   21 P(J)=(Y(J)-SP*V(J))/H
  215 KSKIP=JSKIP
      I=IP1
      IF(I.LE.NM1) GO TO 4
      A(N)=CKJ
      B(NM1)=-B(NM1)
      B(N)=0.D0
      U=DABS(A(1))+DABS(B(1))
      DO 22 I=2,N
   22 U=DMAX1(U,DABS(A(I))+DABS(B(I))+DABS(B(I-1)))
      BD=U
      RBD=1.D0/U
      DO 23 I=1,N
      W(I)=B(I)
      B(I)=(B(I)/U)**2
      A(I)=A(I)/U
      V(I)=0.D0
   23 EV(I)=-1.D0
      U=1.D0
      IK=1
      NDIM=KJ
 1000 K=IK
      EL=EV(K)
   24 ELAM=.5D0*(U+EL)
      DU=(4.D0*DABS(ELAM)+RBD)*EPS
      IF(DABS(U-EL).LE.DU) GO TO 42
      IAG=0
      Q=A(1)-ELAM
      IF(Q.GE.0.D0) IAG=IAG+1
      DO 38 I=2,N
      IF(Q.EQ.0.D0) X=DABS(W(I-1)/BD)/EPS
      IF(Q.NE.0.D0) X=B(I-1)/Q
      Q=A(I)-ELAM-X
      IF( Q.GE.0.D0) IAG=IAG+1
   38 CONTINUE
      IF(IAG.GE.K) GO TO 39
      U=ELAM
      GO TO 24
   39 IF(IAG.EQ.K) GO TO 41
      M=K+1
      DO 40 MM=M,IAG
   40 EV(MM)=ELAM
   41 EL=ELAM
      GO TO 24
   42 ELAM=BD*ELAM
      EV(K)=ELAM
      IF(IK.EQ.1) GO TO 44
      IF(ELAM.GE.EV(IK-1)) EV(IK)=UMEPS*EV(IK-1)
   44 I=IK
      II=1
      L=N-1
      DO 49 J=1,N
   49 Y(J)=1.D0
   50 DO 51 K=1,N
      P(K)=0.D0
      TB(K)=W(K)
   51 TA(K)=BD*A(K)-EV(I)
      J=1
      DO 57 JP1=2,N
      IF(DABS(TA(J)).LT.DABS(W(J))) GO TO 53
      IF(TA(J).EQ.0.D0) TA(J)=EPS
      F=W(J)/TA(J)
      GO TO 55
   53 F=TA(J)/W(J)
      TA(J)=W(J)
      T=TA(JP1)
      TA(JP1)=TB(J)
      TB(J)=T
      P(J)=TB(JP1)
      TB(JP1)=0.D0
      T=Y(J)
      Y(J)=Y(JP1)
      Y(JP1)=T
   55 TB(JP1)=TB(JP1)-F*P(J)
      TA(JP1)=TA(JP1)-F*TB(J)
      Y(JP1)=Y(JP1)-F*Y(J)
   57 J=JP1
      IF(TA(N).EQ.0.D0) TA(N)=EPS
      IF(TA(L).EQ.0.D0) TA(L)=EPS
      Y(N)=Y(N)/TA(N)
      Y(L)=(Y(L)-Y(N)*TB(L))/TA(L)
      DO 62 J=2,L
      K=N-J
      IF(TA(K).EQ.0.D0) TA(K)=EPS
   62 Y(K)=(Y(K)-Y(K+1)*TB(K)-Y(K+2)*P(K))/TA(K)
      AY=DABS(Y(1))
      DO 63 J=2,N
   63 AY=DMAX1(AY,DABS(Y(J)))
      DO 64 J=1,N
   64 Y(J)=Y(J)/AY
      II=II+1
      IF(II.LE.2) GO TO 50
      ID=NDIM-2
      L=N-2
      DO 68 J=1,L
      ID=ID-J-2
      M=N-J
      H=W(M-1)
      IF(H.EQ.0.D0) GO TO 68
      JP1=J+1
      T=DOT(C(ID+1),1,Y(M),1,JP1)/(H*C(ID+1))
      KJ=ID
      DO 67 K=M,N
      KJ=KJ+1
   67 Y(K)=Y(K)+T*C(KJ)
   68 CONTINUE
      XNORM=DSQRT(DOT(Y,1,Y,1,N))
      DO 70 J=1,N
   70 Y(J)=Y(J)/XNORM
C     WRITE(IU) EV(IK),(Y(J),J=1,N)
C     WRITE(6,901) IK,EV(IK)
C 901 FORMAT(I6,5X,D20.12)
      DO 90 J=1,N
   90 UU(J,IK)=Y(J)
      IK=IK+1
      IF(IK.LE.N) GO TO 1000
C     END FILE IU
C     REWIND IU
      RETURN
      END
