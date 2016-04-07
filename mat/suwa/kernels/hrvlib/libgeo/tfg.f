C     *** TFG ********* VERSION 1, MODIFICATION LEVEL 0 *** DKO20812 ***
C     *                                                                *
C     *   FAST FOURIER TRANSFORM , MULTIPLE RADIX FORM                 *
C     *                                                                *
C     *   5736-XM7 COPYRIGHT IBM CORP. 1971                            *
C     *   REFER TO INSTRUCTIONS ON COPYRIGHT NOTICE FORM NO. 120-2083  *
C     *   FE SERVICE NO. 200281                                        *
C     *                                                                *
C     ******************************************************************
C
      SUBROUTINE TFG(RA,CA,NA,IOPT,IER,AUX,NMAX)
      DIMENSION RA(1),CA(1),AUX(2),NF(15)
      DOUBLE PRECISION CD,CE,DA,DC,DN,DV,RD,RT
      EQUIVALENCE (HX,IX)
c     equivalence (DN,SN)
C
C     CHECK FOR TERMINATING ERRORS
C
      ISW  =IER
      IER  =0
      IF (NA-1) 1,1,2
    1 IER  =1000
      GOTO 109
    2 IF (IOPT) 4,3,4
    3 IER  =2000
      GOTO 109
    4 NN   =NA
C
C     DETERMINE THE FACTORS OF NA
C
      M    =0
      JF   =0
      KS   =0
      KA   =4
      KC   =4
      KB   =16
      KD   =-2
    5 IF (NN-KC) 14,6,6
    6 KK   =NN/KB
      IF (NN-KK*KB) 7,13,7
    7 KA   =KA+KD
      IF (KD) 8,9,9
    8 KD   =1
      GOTO 10
    9 KD   =2
   10 KB   =KA
      IF (KS) 11,11,12
   11 KB   =KB*KB
   12 KC   =KB
      GOTO 5
   13 M    =M+1
      JF   =KA
      NF(M)=KA
      NN   =KK
      GOTO 5
   14 IF (KS) 15,15,16
   15 MS   =M
      KA   =2
      KB   =2
      KC   =2
      KD   =1
      KS   =1
      MF   =JF
      MP   =NN-1
      GOTO 5
   16 L    =MS
      IF (M-L-1) 17,17,18
   17 MP   =M+M
   18 IF (L) 20,20,19
   19 M    =M+1
      NF(M)=NF(L)
      L    =L-1
      GOTO 18
C
C     DETERMINE AUXILIARY STORAGE SIZE
C
   20 IF (MF-JF) 21,22,22
   21 MF   =JF
   22 IF (MF-5) 23,23,24
   23 MF   =1
   24 IF (MP-M) 25,25,26
   25 MP   =M+1
   26 K    =MP+4*MF
      IF (NMAX-K) 27,28,28
   27 NMAX =K
      IER  =3000
      GOTO 109
   28 NMAX =K
C
C     INITIALIZE FOURIER TRANSFORM
C
      NJ   =IABS(IOPT)
      RD   =6.2831853071795865D0
      SF   =1.2566370614359173D0
      CF   =COS(SF)
      SF   =SIN(SF)
      SS   =8.660254037844386D-1
      IF (IOPT) 29,30,30
   29 SF   =-SF
      SS   =-SS
      RD   =-RD
   30 CV   =CF*CF-SF*SF
      SV   =CF*(SF+SF)
      RT   =NJ*RD
      NT   =NJ*NA
      KS   =NT
      NN   =KS-NJ
      I    =0
      LK   =0
C
C     START TRANSFORM FOR NEXT FACTOR
C
   31 DA   =RT/KS
      CD   =DCOS(DA)
      SD   =DSIN(DA)
      DA   =0.D0
      DV   =0.D0
      DC   =1.D0
      I    =I+1
      K    =NF(I)
      KL   =KS
      KS   =KS/K
      KK   =1
      KA   =KS+2
      GOTO (28,32,38,40,49),K
      GOTO 51
C
C     PERFORM FOURIER TRANSFORM FOR FACTOR 2
C     (INCL. MULTIPLICATION BY ROTATION FACTOR)
C
   32 KB   =KK+KS
      AJ   =RA(KK)
      AK   =RA(KB)
      RA(KK)=AJ+AK
      RA(KB)=AJ-AK
      BJ   =CA(KK)
      BK   =CA(KB)
      CA(KK)=BJ+BK
      CA(KB)=BJ-BK
      KK   =KB+KS
      IF (KK-NN) 32,32,33
   33 KK   =KK-NN
      IF (KK-KS) 34,34,70
   34 CZ   =CD
      SA   =SD
      DC   =CD
      DV   =CD
      DA   =1.D0
   35 KB   =KK+KS
      AJ   =RA(KK)
      AK   =RA(KB)
      BJ   =CA(KK)
      BK   =CA(KB)
      RA(KK)=AJ+AK
      CA(KK)=BJ+BK
      AK   =AJ-AK
      BK   =BJ-BK
      RA(KB)=CZ*AK-SA*BK
      CA(KB)=SA*AK+CZ*BK
      KK   =KB+KS
      IF (KK-NN) 35,35,36
   36 KB   =KK-NT
      CZ   =-CZ
      KK   =KA-KB
      IF (KK-KB) 37,37,35
   37 DN   =DV+DC
      sn=dn
      SA   =SN*SD
      DV   =DN*CD
      DC   =DV-DA
      CZ   =DC
      DA   =DN
      KK   =KK+NJ
      IF (KK-KB) 35,31,31
C
C     PERFORM FOURIER TRANSFORM FOR FACTOR 3
C
   38 KA   =KK+KS
      KB   =KA+KS
      AK   =RA(KK)
      AL   =RA(KA)
      AM   =RA(KB)
      BK   =CA(KK)
      BL   =CA(KA)
      BM   =CA(KB)
      AJ   =AL+AM
      BJ   =BL+BM
      RA(KK)=AK+AJ
      CA(KK)=BK+BJ
      AK   =AK-0.5E0*AJ
      BK   =BK-0.5E0*BJ
      AJ   =SS*(AL-AM)
      BJ   =SS*(BL-BM)
      RA(KA)=AK-BJ
      RA(KB)=AK+BJ
      CA(KA)=BK+AJ
      CA(KB)=BK-AJ
      KK   =KB+KS
      IF (KK-NN) 38,39,39
   39 KK   =KK-NN
      IF (KK-KS) 38,38,64
C
C     PERFORM FOURIER TRANSFORM FOR FACTOR 4
C     (INCL. MULTIPLICATION BY ROTATION FACTOR)
C
   40 KA   =KK+KS
      KB   =KA+KS
      KC   =KB+KS
      AQ   =RA(KK)
      AP   =RA(KA)
      AM   =RA(KB)
      AL   =RA(KC)
      AK   =AQ+AM
      AM   =AQ-AM
      AJ   =AP+AL
      AL   =AP-AL
      RA(KK)=AK+AJ
      AJ   =AK-AJ
      BQ   =CA(KK)
      BP   =CA(KA)
      BM   =CA(KB)
      BL   =CA(KC)
      BK   =BQ+BM
      BM   =BQ-BM
      BJ   =BP+BL
      BL   =BP-BL
      CA(KK)=BK+BJ
      BJ   =BK-BJ
      IF (IOPT) 41,41,42
   41 AK   =AM+BL
      AM   =AM-BL
      BK   =BM-AL
      BM   =BM+AL
      GOTO 43
   42 AK   =AM-BL
      AM   =AM+BL
      BK   =BM+AL
      BM   =BM-AL
   43 IF (DA) 44,45,44
   44 RA(KA)=CZ*AK-SA*BK
      RA(KB)=CB*AJ-SB*BJ
      RA(KC)=CC*AM-SC*BM
      CA(KA)=SA*AK+CZ*BK
      CA(KB)=SB*AJ+CB*BJ
      CA(KC)=SC*AM+CC*BM
      GOTO 46
   45 RA(KA)=AK
      RA(KB)=AJ
      RA(KC)=AM
      CA(KA)=BK
      CA(KB)=BJ
      CA(KC)=BM
   46 KK   =KC+KS
      IF (KK-NN) 40,47,47
   47 DN   =DV+DC
      sn=dn
      SA   =SN*SD
      DV   =DN*CD
      DC   =DV-DA
      CZ   =DC
      DA   =DN
      CC   =CZ+CZ
      SB   =CC*SA
      SC   =CC*CZ
      CB   =-0.5E0+SC-0.5E0
      SC   =SC+CB
      CC   =SC*CZ-CC
      SC   =SC*SA
      KK   =KK-NN
      IF (KK-KS) 40,40,48
   48 IF (KS-NJ) 70,70,31
C
C     PERFORM FOURIER TRANSFORM FOR FACTOR 5
C
   49 KA   =KK+KS
      KB   =KA+KS
      KC   =KB+KS
      KD   =KC+KS
      AN   =RA(KA)
      AM   =RA(KD)
      AP   =AN+AM
      AM   =AN-AM
      AN   =RA(KB)
      AL   =RA(KC)
      AQ   =AN+AL
      AL   =AN-AL
      AN   =RA(KK)
      RA(KK)=AN+AP+AQ
      BN   =CA(KA)
      BM   =CA(KD)
      BP   =BN+BM
      BM   =BN-BM
      BN   =CA(KB)
      BL   =CA(KC)
      BQ   =BN+BL
      BL   =BN-BL
      BN   =CA(KK)
      CA(KK)=BN+BP+BQ
      AK   =CF*AP+CV*AQ+AN
      BK   =CF*BP+CV*BQ+BN
      AJ   =SF*AM+SV*AL
      BJ   =SF*BM+SV*BL
      RA(KA)=AK-BJ
      RA(KD)=AK+BJ
      CA(KA)=BK+AJ
      CA(KD)=BK-AJ
      AK   =CV*AP+CF*AQ+AN
      BK   =CV*BP+CF*BQ+BN
      AJ   =SV*AM-SF*AL
      BJ   =SV*BM-SF*BL
      RA(KB)=AK-BJ
      RA(KC)=AK+BJ
      CA(KB)=BK+AJ
      CA(KC)=BK-AJ
      KK   =KD+KS
      IF (KK-NN) 49,50,50
   50 KK   =KK-NN
      IF (KK-KS) 49,49,64
C
C     PERFORM FOURIER TRANSFORM FOR FACTORS 7,11,13,...
C
   51 IF (LK-K) 52,54,52
   52 LK   =K
      JB   =MP+MF
      JC   =JB+MF
      JS   =JC+MF
      JD   =JC+K
      JT   =JS+K
      DN   =RD/K
      sn=dn
      CE   =DCOS(DN)
      SA   =DSIN(DN)
      AUX(JD)=1.
      AUX(JT)=0.
      KA   =(K-1)/2
      DO 53 L=1,KA
         DN   =DV+DC
         sn=dn
         SB   =SN*SA
         DV   =DN*CE
         DC   =DV-DA
         DA   =DN
         J    =JC+L
         AUX(J)=DC
         J    =JS+L
         AUX(J)=SB
         J    =JD-L
         AUX(J)=DC
         J    =JT-L
   53    AUX(J)=-SB
   54 AN   =RA(KK)
      BN   =CA(KK)
      AK   =AN
      BK   =BN
      KA   =KK+KS
      KB   =KK+KL
      L    =MP+1
      J    =JB+1
   55 KB   =KB-KS
      AQ   =RA(KA)
      AP   =RA(KB)
      BQ   =CA(KA)
      BP   =CA(KB)
      AL   =AQ+AP
      BL   =BQ+BP
      L    =L+1
      J    =J+1
      AUX(L)=AL
      AUX(J)=BL
      AK   =AL+AK
      BK   =BL+BK
      L    =L+1
      J    =J+1
      AUX(L)=AQ-AP
      AUX(J)=BQ-BP
      KA   =KA+KS
      IF (KA-KB) 55,56,56
   56 RA(KK)=AK
      CA(KK)=BK
      KA   =KK
      KB   =KK+KL
      KD   =1
   57 KA   =KA+KS
      KB   =KB-KS
      JK   =KD
      AK   =AN
      BK   =BN
      AJ   =0.E0
      BJ   =0.E0
      JL   =1
   58 JL   =JL+1
      K    =JC+JK
      L    =MP+JL
      J    =JB+JL
      AK   =AUX(K)*AUX(L)+AK
      BK   =AUX(K)*AUX(J)+BK
      JL   =JL+1
      K    =JS+JK
      L    =L+1
      J    =J+1
      AJ   =AUX(K)*AUX(L)+AJ
      BJ   =AUX(K)*AUX(J)+BJ
      JK   =JK+KD
      IF (JK-LK) 60,60,59
   59 JK   =JK-LK
   60 IF (JL-LK) 58,61,61
   61 JL   =LK-KD
      RA(KA)=AK-BJ
      RA(KB)=AK+BJ
      CA(KA)=BK+AJ
      CA(KB)=BK-AJ
      KD   =KD+1
      IF (KD-JL) 57,62,62
   62 KK   =KK+KL
      IF (KK-NN) 54,54,63
   63 KK   =KK-NN
      IF (KK-KS) 54,54,64
C
C     MULTIPLY ODD-FACTOR TRANSFORM BY ROTATION FACTOR
C
   64 IF (I-M) 65,70,70
   65 KK   =NJ+1
      KA   =NT-KS
      KB   =KL-NJ
      DA   =1.D0
      DC   =CD
      DV   =CD
      CB   =CD
      SA   =SD
   66 CZ   =CB
      CC   =CB
      SB   =SA
      SC   =1.E0
      KK   =KK+KS
   67 AK   =RA(KK)
      BK   =CA(KK)
      RA(KK)=CB*AK-SB*BK
      CA(KK)=SB*AK+CB*BK
      KK   =KK+KL
      IF (KK-NT) 67,67,68
   68 SN   =CC+CB
      SB   =SN*SA
      CC   =SN*CZ
      CB   =CC-SC
      SC   =SN
      KK   =KK-KA
      IF (KK-KL) 67,67,69
   69 DN   =DV+DC
      sn=dn
      SA   =SN*SD
      DV   =DN*CD
      DC   =DV-DA
      CB   =DC
      DA   =DN
      KK   =KK-KB
      IF (KK-KS) 66,66,31
C
C     PERMUTE VECTORS RA,CA TO ORIGINAL ORDER
C
   70 KD   =MS+MS+1
      K    =KD+2
      IF (M-KD) 71,72,72
   71 K    =K-1
   72 KA   =NT
      KB   =1
      KC   =NJ
      KK   =MS+1
      DO 73 L=1,KK
         KA   =KA/KB
         KC   =KC*KB
         KB   =NF(L)
         IX   =KA
         AUX(L)=HX
         J    =K-L
         IX   =KC
   73    AUX(J)=HX
      KL   =KA
      JF   =KC
      IF (MS) 81,81,74
C
C     PERMUTATIONS CORRESPONDING TO SQUARE FACTORS OF NA
C
   74 HX   =AUX(2)
      KS   =IX
      KK   =NJ+1
      KB   =KS+1
   75 IF (KK-KB) 76,77,77
   76 AK   =RA(KK)
      RA(KK)=RA(KB)
      RA(KB)=AK
      BK   =CA(KK)
      CA(KK)=CA(KB)
      CA(KB)=BK
   77 KK   =KK+NJ
      KB   =KB+KS
      IF (KB-NT) 75,78,78
   78 IF (KK-NT) 79,81,81
   79 KA   =NT
      KC   =KS
      L    =2
   80 L    =L+1
      KB   =KB-KA
      KA   =KC
      HX   =AUX(L)
      KC   =IX
      KB   =KB+KC
      IF (KB-KA) 75,75,80
C
C     PERMUTATIONS CORRESP. TO SQUARE-FREE FACTORS OF NA
C
   81 IF (KD-M) 82,111,111
   82 KA   =1
      L    =M-MS+1
      GOTO 84
   83 KD   =KA
      KA   =KA*NF(L)
   84 NF(L)=KA
      L    =L-1
      IF (L-MS) 85,85,83
   85 NN   =KA-1
      MS   =MS+2
      KC   =0
      DO 89 L=1,NN
         K    =MS
         KB   =KA
         KK   =KD
   86    KC   =KC+KK
         IF (KC-KB) 88,87,87
   87    KC   =KC-KB
         KB   =KK
         K    =K+1
         KK   =NF(K)
         GOTO 86
   88    IX   =KC
   89    AUX(L)=HX
      DO 94 L=1,NN
         HX   =AUX(L)
         KK   =IX
         IF (KK) 94,94,90
   90    IF (KK-L) 91,93,91
   91    K    =KK
         HX   =AUX(K)
         KK   =IX
         IX   =-IX
         AUX(K)=HX
         IF (KK-L) 91,92,91
   92    KC   =KK
         GOTO 94
   93    IX   =-L
         AUX(L)=HX
   94    CONTINUE
      MF   =NJ*MF
      GOTO 108
   95 L    =KC
   96 IF (L) 108,108,97
   97 HX   =AUX(L)
      IF (IX) 107,107,98
   98 KD   =JF
   99 KS   =MIN0(KD,MF)
      KD   =KD-KS
      HX   =AUX(L)
      K    =IX
      KK   =K*JF+I+KD
      KA   =KK+KS
      KB   =MP+1
  100 AUX(KB)=RA(KA)
      KB   =KB+1
      AUX(KB)=CA(KA)
      KA   =KA-NJ
      KB   =KB+1
      IF (KA-KK) 101,101,100
  101 HX   =AUX(K)
      NN   =IX
      KA   =KK+KS
      KB   =KA-JF*(K+NN)
      K    =-NN
  102 RA(KA)=RA(KB)
      CA(KA)=CA(KB)
      KA   =KA-NJ
      KB   =KB-NJ
      IF (KA-KK) 103,103,102
  103 KK   =KB
      IF (K-L) 101,104,101
  104 KA   =KK+KS
      KB   =MP+1
  105 RA(KA)=AUX(KB)
      KB   =KB+1
      CA(KA)=AUX(KB)
      KA   =KA-NJ
      KB   =KB+1
      IF (KA-KK) 106,106,105
  106 IF (KD) 99,107,99
  107 L    =L-1
      GOTO 96
  108 NT   =NT-KL
      I    =NT-NJ+1
      IF (NT) 111,95,95
  109 IF (ISW+12345) 110,111,110
  110 CALL WIER(IER,20812)
  111 RETURN
      END
